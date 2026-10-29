"""The synthetic fixture through the processor (both input layouts), the method's pieces, and the Lambda handler."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import wasserstein_distance

from processor import assess as A
from processor import handler
from processor import main as app
from processor.inputs import read_prep_missing
from tests.make_fixture import make_fixture

OUTPUTS = {"emd_reference.csv", "emd_batches.csv", "flags.csv", "summary.json", "provenance.json",
           "heatmap_reference.png", "heatmap_batches.png"}


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, value)
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_shift_is_flagged(fixture_root, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixture_root[layout], out, monkeypatch, REFERENCE_COLUMN="control", REFERENCE_VALUES="ref") == 0
    # Everything under OUTPUT_DIR/<LABEL>/, nothing beside it.
    assert [p.name for p in out.iterdir()] == ["raw"]
    files = {p.name for p in (out / "raw").iterdir()}
    assert OUTPUTS <= files
    assert any(f.startswith("density_M1") for f in files)

    ref = rows(out / "raw" / "emd_reference.csv")
    assert list(ref[0])[:3] == ["batch", "marker", "emd"]
    assert len(ref) == 4 * 5   # 4 batches x 5 markers, one reference value
    flagged_ref = {(r["batch"], r["marker"]) for r in ref if r["flagged"] == "true"}
    assert flagged_ref == {("b3", "M1"), ("b3", "M3")}

    bat = rows(out / "raw" / "emd_batches.csv")
    assert list(bat[0])[:4] == ["batch", "condition", "marker", "emd"]
    assert len(bat) == 4 * 2 * 5
    flagged_bat = {(r["batch"], r["condition"], r["marker"]) for r in bat if r["flagged"] == "true"}
    assert flagged_bat == {("b3", c, m) for c in ("stim", "unstim") for m in ("M1", "M3")}
    # Stimulation raises M2 in every batch alike: biology, so no batch is flagged for it.
    assert not any(r["marker"] == "M2" and r["flagged"] == "true" for r in bat)

    flags = rows(out / "raw" / "flags.csv")
    missing = [f for f in flags if f["comparison"] == "missing"]
    assert [(f["batch"], f["marker"]) for f in missing] == [("b4", "M6")]
    assert {f["batch"] for f in flags if f["level"] == "batch"} == {"b3", "b4"}

    summary = json.loads((out / "raw" / "summary.json").read_text())
    assert summary["label"] == "raw" and summary["n_batches"] == 4 and summary["n_markers"] == 5
    assert summary["median_emd_source"] == "reference"
    scores = {b: e["score"] for b, e in summary["per_batch"].items()}
    assert max(scores, key=scores.get) == "b3"
    assert summary["per_batch"]["b3"]["flagged_markers"] == ["M1", "M3"]

    prov = json.loads((out / "raw" / "provenance.json").read_text())
    assert prov["app"] == "processor-cytometry-batch-assess" and prov["version"] == "0.1.0"
    assert prov["n_ok"] == 20 and prov["n_failed"] == 0
    assert {"numpy", "scipy", "flowio"} <= set(prov["tools"])
    assert all(s["status"] == "ok" for s in prov["samples"])


def test_labels_share_an_output_folder(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, LABEL="raw") == 0
    assert run(fixture_root["linked"], out, monkeypatch, LABEL="cytonorm") == 0
    assert sorted(p.name for p in out.iterdir()) == ["cytonorm", "raw"]
    # Without reference samples: no reference rows, the headline comes from the batch comparison.
    assert rows(out / "raw" / "emd_reference.csv") == []
    summary = json.loads((out / "cytonorm" / "summary.json").read_text())
    assert summary["label"] == "cytonorm" and summary["median_emd_source"] == "batches"
    assert not (out / "raw" / "heatmap_reference.png").exists()


@pytest.mark.parametrize("label", ["Raw", "raw_1", "a/b", "raw.1"])
def test_bad_label(fixture_root, tmp_path, monkeypatch, label):
    assert run(fixture_root["linked"], tmp_path / "out", monkeypatch, LABEL=label) == 2


def test_reference_parameters_go_together(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root["linked"], tmp_path / "out", monkeypatch, REFERENCE_COLUMN="control") == 2


def test_a_bad_file_is_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    (src / "fcs" / "b1_stim1.fcs").write_bytes(b"not an fcs file")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    prov = json.loads((out / "raw" / "provenance.json").read_text())
    assert prov["n_failed"] == 1 and prov["n_ok"] == 19


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_subsample_is_seeded_per_sample():
    a = A.subsample_rows(1000, 100, np.random.default_rng(A.sample_seed(42, "s1")))
    b = A.subsample_rows(1000, 100, np.random.default_rng(A.sample_seed(42, "s1")))
    c = A.subsample_rows(1000, 100, np.random.default_rng(A.sample_seed(42, "s2")))
    assert np.array_equal(a, b) and not np.array_equal(a, c)
    assert A.subsample_rows(50, 100, np.random.default_rng(0)) is None


def test_histogram_emd_matches_exact_emd():
    rng = np.random.default_rng(1)
    x, y = rng.normal(1.0, 0.5, 20000), rng.normal(1.3, 0.7, 15000)
    exact = wasserstein_distance(x, y)
    binned = A.emd(A.normalised(A.histogram(x)[0]), A.normalised(A.histogram(y)[0]))
    assert abs(binned - exact) <= A.BIN_WIDTH


def test_robust_threshold():
    med, mad, thr = A.robust_threshold([0.1, 0.2, 0.3, 0.4, 5.0], 3)
    assert med == pytest.approx(0.3)
    assert mad == pytest.approx(1.4826 * 0.1)
    assert thr == pytest.approx(0.3 + 3 * 1.4826 * 0.1)


def test_prep_json_missing(tmp_path):
    prep = {"samples": {
        "a1": {"batch": "a", "channels_dropped": [{"channel": "X-A", "marker": "GITR"}, {"channel": "Y-A", "marker": "CD28"}]},
        "b1": {"batch": "b", "channels_dropped": [{"channel": "Z-A", "marker": "CD28"}]},
        "c1": {"batch": "c", "channels_dropped": [{"channel": "Y-A", "marker": "CD28"}]},
    }}
    path = tmp_path / "prep.json"
    path.write_text(json.dumps(prep))
    # GITR was only in a1; CD28 is in every sample (on another channel in b1): not missing.
    assert read_prep_missing(path, {}) == {"GITR": {"b1", "c1"}}


def test_handler(fixture_root, tmp_path, monkeypatch):
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    out = tmp_path / "out"
    result = handler.handler({"inputDir": str(fixture_root["linked"]), "outputDir": str(out), "LABEL": "cycombine",
                              "REFERENCE_COLUMN": "control", "REFERENCE_VALUES": ["ref"], "CELLS_PER_SAMPLE": 1000}, None)
    assert result["status"] == "ok"
    summary = json.loads((out / "cycombine" / "summary.json").read_text())
    assert summary["reference"]["n_samples"] == 4
    assert summary["parameters"]["CELLS_PER_SAMPLE"] == 1000
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(fixture_root["linked"]), "outputDir": str(out), "LABEL": "BAD"}, None)
