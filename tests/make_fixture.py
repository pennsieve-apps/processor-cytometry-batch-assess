"""
Fixtures that follow the workflow contract (fcs/<sample_id>.fcs with arcsinh-
transformed marker channels, metadata/{samples.csv, panel.csv, prep.json}).

  python -m tests.make_fixture data/fixture                 synthetic cohort (linked + flat layouts)
  python -m tests.make_fixture data/fixture-demo --from DIR [--events 20000] [--shift batch03:CD3,CD4,CD8:0.4]
      subsample a prep-shaped folder (e.g. workflow-immune-cytometry/data/work/prep) and plant a shift

The synthetic cohort: 4 batches x {unstim, stim} x 2 samples + one reference
sample (control=ref) per batch, 6 markers (bimodal), stim raises M2 (biology),
batch b3 is shifted by +0.6 on M1 and M3 (the planted batch effect), and M6 (on
channel V6) was dropped by prep because batch b4 lacked it (prep.json).
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np

MARKERS = [("B1-A", "M1"), ("B2-A", "M2"), ("B3-A", "M3"), ("B4-A", "M4"), ("B5-A", "M5")]
SCATTER = ["FSC-A", "SSC-A"]
BATCHES = ["b1", "b2", "b3", "b4"]
SHIFT_BATCH, SHIFT_MARKERS, SHIFT = "b3", ("M1", "M3"), 0.6


def write_fcs(path: Path, pnn: list[str], pns: list[str], events: np.ndarray, date: str = "01-Jan-2024") -> None:
    import flowio

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        flowio.create_fcs(fh, events.astype(np.float32).ravel().tolist(), pnn, opt_channel_names=pns,
                          metadata_dict={"DATE": date})


def _synthetic_events(rng: np.random.Generator, n: int, batch: str, condition: str) -> np.ndarray:
    cols = [rng.normal(1e5, 2e4, n), rng.normal(5e4, 1e4, n)]
    for i, (_, marker) in enumerate(MARKERS):
        pos_frac = 0.3 + 0.08 * i + (0.3 if (marker == "M2" and condition == "stim") else 0.0)
        pos = rng.random(n) < pos_frac
        values = np.where(pos, rng.normal(2.5 + 0.2 * i, 0.45, n), rng.normal(0.0, 0.35, n))
        values += rng.normal(0, 0.03)   # small, real batch-to-batch noise
        if batch == SHIFT_BATCH and marker in SHIFT_MARKERS:
            values += SHIFT
        cols.append(values)
    cols.append(np.sort(rng.uniform(0, 1000, n)))   # Time
    return np.column_stack(cols)


def make_synthetic(root: Path, events: int = 3000, seed: int = 7) -> Path:
    rng = np.random.default_rng(seed)
    pnn = SCATTER + [c for c, _ in MARKERS] + ["Time"]
    pns = ["", ""] + [m for _, m in MARKERS] + [""]
    rows, prep_samples = [], {}
    for b in BATCHES:
        plan = [("ref", "unstim", "ref")] + [(f"{c}{k}", c, "none") for c in ("unstim", "stim") for k in (1, 2)]
        for name, condition, control in plan:
            sid = f"{b}_{name}"
            write_fcs(root / "fcs" / f"{sid}.fcs", pnn, pns, _synthetic_events(rng, events, b, condition))
            rows.append({"sample_id": sid, "file": f"{sid}_raw.fcs", "batch": b, "condition": condition,
                         "n_events": events, "control": control})
            had_v6 = b != "b4"
            prep_samples[sid] = {"batch": b, "channels_dropped": (
                [{"channel": "V6-A", "marker": "M6", "reason": "missing in 5 of 20 samples"}] if had_v6 else [])}
    meta = root / "metadata"
    meta.mkdir(parents=True, exist_ok=True)
    with open(meta / "samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(meta / "panel.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["channel", "marker", "type", "cofactor", "clustering", "state"])
        for c in SCATTER:
            w.writerow([c, c, "scatter", "", "false", "false"])
        for c, m in MARKERS:
            w.writerow([c, m, "marker", "6000", "true", "false"])
        w.writerow(["Time", "Time", "time", "", "false", "false"])
    (meta / "prep.json").write_text(json.dumps({"app": "processor-fcs-prep", "channelsDropped": ["V6-A"],
                                                "samples": prep_samples}, indent=1))
    return root


def flatten(linked: Path, flat: Path) -> Path:
    """The same files with no port folders, nested a little, as an unlinked run would see them."""
    for p in linked.rglob("*"):
        if p.is_file():
            dest = flat / ("nested" if p.suffix == ".fcs" else "") / p.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(p, dest)
    return flat


def make_fixture(root: Path, layout: str = "both") -> dict[str, Path]:
    root = Path(root)
    out = {"linked": make_synthetic(root / "linked")}
    if layout == "both":
        out["flat"] = flatten(out["linked"], root / "flat")
    return out


def subsample_prep(source: Path, root: Path, events: int, shift: str | None, seed: int = 42) -> Path:
    """Subsample a prep-shaped folder (fcs/ + metadata/) and plant a shift: 'batch:M1,M2:delta'."""
    from processor.fcs import read_events, read_header

    shift_batch, shift_markers, delta = None, set(), 0.0
    if shift:
        shift_batch, markers, d = shift.split(":")
        shift_markers, delta = set(markers.split(",")), float(d)
    meta_src = source / "metadata"
    (root / "metadata").mkdir(parents=True, exist_ok=True)
    for name in ("samples.csv", "panel.csv", "prep.json"):
        if (meta_src / name).exists():
            shutil.copy(meta_src / name, root / "metadata" / name)
    with open(meta_src / "samples.csv", newline="") as fh:
        batch_of = {r["sample_id"]: r["batch"] for r in csv.DictReader(fh)}
    rng = np.random.default_rng(seed)
    for path in sorted((source / "fcs").glob("*.fcs")):
        header = read_header(path)
        _, ev = read_events(path)
        if ev.shape[0] > events:
            ev = ev[np.sort(rng.choice(ev.shape[0], events, replace=False))]
        if batch_of.get(path.stem) == shift_batch:
            for j, s in enumerate(header.pns):
                if s in shift_markers:
                    ev[:, j] += delta
        write_fcs(root / "fcs" / path.name, header.pnn, header.pns, ev, header.keywords.get("date", "01-Jan-2024"))
        print(f"{path.stem}: {ev.shape[0]} events{' (shifted)' if batch_of.get(path.stem) == shift_batch else ''}")
    return root


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--layout", default="both", choices=["both", "linked"])
    ap.add_argument("--from", dest="source")
    ap.add_argument("--events", type=int, default=20000)
    ap.add_argument("--shift", default=None, help="batch:M1,M2:delta, e.g. batch03:CD3,CD4,CD8:0.4")
    a = ap.parse_args()
    if a.source:
        subsample_prep(Path(a.source), Path(a.root), a.events, a.shift)
    else:
        print(make_fixture(Path(a.root), a.layout))
