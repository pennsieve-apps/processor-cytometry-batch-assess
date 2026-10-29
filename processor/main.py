"""
processor-cytometry-batch-assess: how much do batches differ, per marker?

Reads the transformed FCS files of a cohort (port `fcs`) and its samples.csv /
panel.csv / prep.json (port `metadata`), and writes, all under
OUTPUT_DIR/<LABEL>/ so several runs (raw, after CytoNorm, after cyCombine) can
be linked into one folder of the report step:

  emd_reference.csv   batch x marker EMD of the reference samples (when named)
  emd_batches.csv     batch x condition x marker EMD against the other batches
  flags.csv           flagged batch x marker distances, missing markers, flagged batches
  summary.json        label, counts, median EMD, per-batch score, flags
  heatmap_*.png, density_<marker>.png, provenance.json

See processor/assess.py for the method.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import platform
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-cytometry-batch-assess")
log.setLevel(logging.INFO)   # the Lambda runtime sets the root logger to WARNING first
for noisy in ("matplotlib", "PIL", "fontTools"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-cytometry-batch-assess"
REPO_ROOT = Path(__file__).resolve().parent.parent
LABEL_RE = re.compile(r"^[a-z0-9-]+$")

DEFAULTS = {
    "LABEL": "raw",
    "REFERENCE_COLUMN": "",
    "REFERENCE_VALUES": "",
    "CONDITION_COLUMN": "condition",
    "CELLS_PER_SAMPLE": "50000",
    "MAD_K": "3",
    "MIN_EMD": "0.1",
    "BATCH_FLAG_FRACTION": "0.2",
    "MARKERS": "",
    "MAX_PLOTS": "6",
    "SEED": "42",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {}
    p["LABEL"] = raw["LABEL"].strip()
    if not LABEL_RE.match(p["LABEL"]):
        raise ParameterError(f"LABEL must match [a-z0-9-]+, got {p['LABEL']!r}")
    p["REFERENCE_COLUMN"] = raw["REFERENCE_COLUMN"].strip()
    p["REFERENCE_VALUES"] = [v.strip() for v in raw["REFERENCE_VALUES"].split(",") if v.strip()]
    if bool(p["REFERENCE_COLUMN"]) != bool(p["REFERENCE_VALUES"]):
        raise ParameterError("REFERENCE_COLUMN and REFERENCE_VALUES go together: set both or neither")
    p["CONDITION_COLUMN"] = raw["CONDITION_COLUMN"].strip()
    p["MARKERS"] = [v.strip() for v in raw["MARKERS"].split(",") if v.strip()]
    try:
        p["CELLS_PER_SAMPLE"] = int(raw["CELLS_PER_SAMPLE"])
        p["MAX_PLOTS"] = int(raw["MAX_PLOTS"])
        p["SEED"] = int(raw["SEED"])
        p["MAD_K"] = float(raw["MAD_K"])
        p["MIN_EMD"] = float(raw["MIN_EMD"])
        p["BATCH_FLAG_FRACTION"] = float(raw["BATCH_FLAG_FRACTION"])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["CELLS_PER_SAMPLE"] < 0 or p["MAD_K"] < 0 or p["MIN_EMD"] < 0 or not 0 <= p["BATCH_FLAG_FRACTION"] <= 1:
        raise ParameterError("CELLS_PER_SAMPLE, MAD_K, MIN_EMD must be >= 0 and BATCH_FLAG_FRACTION in [0, 1]")
    return p


def _rel(path: Path | None, root: Path) -> str | None:
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _fmt(v):
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.6g}"
    if isinstance(v, bool):
        return "true" if v else "false"
    return v


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r.get(c, "")) for c in columns])


def _clean_json(obj):
    if isinstance(obj, float):
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, dict):
        return {k: _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    import flowio
    import matplotlib
    import numpy
    import pandas
    import scipy

    return {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
            "pandas": pandas.__version__, "matplotlib": matplotlib.__version__, "flowio": flowio.__version__}


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)   # bytes on macOS, KiB on Linux
    except (ImportError, OSError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2

    from processor import assess as A
    from processor import figures as F
    from processor.fcs import read_events
    from processor.inputs import discover

    try:
        found = discover(input_dir, "batch", params["CONDITION_COLUMN"])
    except (OSError, ValueError) as exc:
        log.error("cannot read the metadata: %s", exc)
        return 1
    for name, reason in found.skipped.items():
        log.warning("%s: skipped, %s", name, reason)
    if not found.samples:
        log.error("no FCS samples found under %s", input_dir)
        return 1
    markers = found.markers
    if params["MARKERS"]:
        wanted = set(params["MARKERS"])
        markers = [(c, m) for c, m in markers if m in wanted or c in wanted]
        unknown = wanted - {m for _, m in markers} - {c for c, _ in markers}
        if unknown:
            log.warning("MARKERS not in the panel: %s", ", ".join(sorted(unknown)))
    if not markers:
        log.error("no marker channels (panel.csv type=marker) to assess")
        return 1
    marker_names = [m for _, m in markers]
    if params["REFERENCE_COLUMN"] and found.samples and params["REFERENCE_COLUMN"] not in found.samples[0].row:
        log.error("REFERENCE_COLUMN %r is not a column of samples.csv", params["REFERENCE_COLUMN"])
        return 1
    log.info("%d samples, %d batches, %d markers, label %s", len(found.samples),
             len({s.batch for s in found.samples}), len(markers), params["LABEL"])

    hists, statuses = [], []
    t_read = time.monotonic()
    for s in found.samples:
        t = time.monotonic()
        status = {"sample_id": s.sample_id, "batch": s.batch, "condition": s.condition,
                  "file": _rel(s.path, input_dir)}
        try:
            channels, events = read_events(s.path)
            h = A.sample_histograms(s.sample_id, s.batch, s.condition, channels.pnn, events, markers,
                                    params["CELLS_PER_SAMPLE"], params["SEED"])
            del events
            if not h.counts:
                raise ValueError("none of the panel's marker channels are in this file")
            hists.append(h)
            clipped = {m: round(f, 5) for m, f in h.clipped.items() if f > 0.001}
            status.update(status="ok", n_events=h.n_events, n_cells=h.n_cells, markers_absent=h.absent,
                          values_outside_grid=clipped, seconds=round(time.monotonic() - t, 3))
            if clipped:
                log.warning("%s: values outside [%g, %g] in %s; are the markers transformed?", s.sample_id,
                            A.GRID_LO, A.GRID_HI, ", ".join(clipped))
            log.info("%s: batch %s, condition %s, %d events, %d used, %d markers%s, %.2fs", s.sample_id, s.batch,
                     s.condition, h.n_events, h.n_cells, len(h.counts),
                     f", lacks {', '.join(h.absent)}" if h.absent else "", time.monotonic() - t)
        except Exception as exc:   # one bad file is skipped, not fatal
            status.update(status="failed", error=f"{type(exc).__name__}: {exc}", seconds=round(time.monotonic() - t, 3))
            log.error("%s: failed: %s", s.sample_id, exc)
        statuses.append(status)
    t_read = time.monotonic() - t_read
    rss_read = peak_rss_mb()
    log.info("read %d samples in %.1fs, peak RSS %s MB", len(statuses), t_read, rss_read)
    if not hists:
        log.error("every sample failed")
        return 1

    t_cmp = time.monotonic()
    batches = sorted({h.batch for h in hists})
    if len(batches) < 2:
        log.warning("only one batch (%s): nothing to compare", batches[0])
    ref_of = {}
    if params["REFERENCE_COLUMN"]:
        values = set(params["REFERENCE_VALUES"])
        ref_of = {s.sample_id: s.row.get(params["REFERENCE_COLUMN"], "NA") for s in found.samples
                  if s.row.get(params["REFERENCE_COLUMN"]) in values}
        ref_batches = {s.batch for s in found.samples if s.sample_id in ref_of}
        log.info("reference samples: %d in %d batches (%s in %s)", len(ref_of), len(ref_batches),
                 ",".join(params["REFERENCE_VALUES"]), params["REFERENCE_COLUMN"])
        if not ref_of:
            log.warning("no sample has %s in %s; skipping the reference comparison",
                        params["REFERENCE_COLUMN"], params["REFERENCE_VALUES"])
    ref_rows = A.compare_reference(hists, ref_of, marker_names) if ref_of else []
    batch_rows = A.compare_batches(hists, marker_names)
    A.mark_outliers(ref_rows, ("reference", "marker"), params["MAD_K"], params["MIN_EMD"])
    A.mark_outliers(batch_rows, ("condition", "marker"), params["MAD_K"], params["MIN_EMD"])
    batch_of = {h.sample_id: h.batch for h in hists}
    missing = A.missing_by_batch(hists, found.missing, batch_of)
    flags, summary = A.summarise(params["LABEL"], batches, marker_names, ref_rows, batch_rows, missing,
                                 params["BATCH_FLAG_FRACTION"])
    t_cmp = time.monotonic() - t_cmp

    out = output_dir / params["LABEL"]
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "emd_reference.csv", ref_rows,
              ["batch", "marker", "emd", "reference", "threshold", "flagged", "n_samples", "n_cells"])
    write_csv(out / "emd_batches.csv", batch_rows,
              ["batch", "condition", "marker", "emd", "threshold", "flagged", "n_samples", "n_cells"])
    write_csv(out / "flags.csv", flags,
              ["level", "comparison", "batch", "condition", "reference", "marker", "emd", "threshold", "reason"])
    for b, e in summary["per_batch"].items():
        log.info("batch %s: score %s (%s), %d flagged markers%s", b,
                 "n/a" if e["score"] is None else f"{e['score']:.3f}", e["score_source"], e["n_flagged_markers"],
                 f" ({', '.join(e['flagged_markers'])})" if e["flagged_markers"] else "")

    t_fig = time.monotonic()
    figures = []
    try:
        figures = make_figures(out, params, batches, marker_names, hists, ref_of, ref_rows, batch_rows, missing,
                               summary, A, F)
    except Exception as exc:   # figures are a convenience; the tables stand without them
        log.exception("figures failed: %s", exc)
    t_fig = time.monotonic() - t_fig

    summary.update({
        "n_samples": len(hists),
        "reference": {"column": params["REFERENCE_COLUMN"] or None, "values": params["REFERENCE_VALUES"],
                      "n_samples": len(ref_of)},
        "condition_column": params["CONDITION_COLUMN"],
        "markers": marker_names,
        "parameters": {k: v for k, v in params.items()},
        "figures": [f.name for f in figures],
    })
    write_json(out / "summary.json", summary)
    n_failed = sum(1 for s in statuses if s["status"] != "ok")
    write_json(out / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(),
        "parameters": params,
        "method": {"emd": "scipy.stats.wasserstein_distance on per-marker histograms",
                   "grid": [A.GRID_LO, A.GRID_HI], "bin_width": A.BIN_WIDTH,
                   "threshold": "median + MAD_K * 1.4826 * MAD over batches, and > MIN_EMD",
                   "pooling": "equal weight per sample within a batch, per batch across batches"},
        "inputs": {"fcs": [_rel(p, input_dir) for p in found.files_seen],
                   "samples_csv": _rel(found.samples_csv, input_dir), "panel_csv": _rel(found.panel_csv, input_dir),
                   "prep_json": _rel(found.prep_json, input_dir), "skipped": found.skipped},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"read": round(t_read, 2), "compare": round(t_cmp, 2), "figures": round(t_fig, 2),
                            "total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": {"after_reading": rss_read, "total": peak_rss_mb()},
        "samples": statuses, "n_ok": len(hists), "n_failed": n_failed,
        "outputs": sorted(p.name for p in out.iterdir()) + ["provenance.json"],
    })
    log.info("done: %d samples (%d failed), %d batches, %d flags, median EMD %s, %.1fs -> %s", len(hists), n_failed,
             len(batches), len(flags), "n/a" if summary["median_emd"] is None else f"{summary['median_emd']:.3f}",
             time.monotonic() - t0, out)
    return 0


GLYPH = {"missing": "\u2715", "only_here": "\u25c6"}


def make_figures(out, params, batches, markers, hists, ref_of, ref_rows, batch_rows, missing, summary, A, F):
    import numpy as np

    paths = []
    all_values = [r["emd"] for r in ref_rows + batch_rows]
    vmax = max(0.25, float(np.ceil(max(all_values) * 10) / 10)) if all_values else 1.0
    # Markers some batches lack (e.g. dropped by fcs-prep) get a column too, marked \u2715 where flagged.
    markers = list(markers) + sorted({m for (_, m) in missing} - set(markers))
    col = {m: j for j, m in enumerate(markers)}

    if ref_rows:
        values = np.full((len(batches), len(markers)), np.nan)
        flagged = np.zeros_like(values, dtype=bool)
        acc = defaultdict(list)
        for r in ref_rows:
            acc[(r["batch"], r["marker"])].append(r)
        for (b, m), rs in acc.items():
            i, j = batches.index(b), col[m]
            values[i, j] = float(np.mean([r["emd"] for r in rs]))
            flagged[i, j] = any(r["flagged"] for r in rs)
        marks = {}
        for (b, m), info in missing.items():
            if m in col and b in batches:
                flagged[batches.index(b), col[m]] = True
                marks[(batches.index(b), col[m])] = GLYPH[info["kind"]]
        refs = ", ".join(params["REFERENCE_VALUES"])
        paths.append(F.heatmap([("", batches, markers, values, flagged, marks)],
                               f"Reference samples ({params['REFERENCE_COLUMN']}={refs}) vs pooled reference · {params['LABEL']}",
                               vmax, out / "heatmap_reference.png"))
    if batch_rows:
        panels = []
        for c in sorted({r["condition"] for r in batch_rows}):
            rows_c = [r for r in batch_rows if r["condition"] == c]
            bs = sorted({r["batch"] for r in rows_c})
            values = np.full((len(bs), len(markers)), np.nan)
            flagged = np.zeros_like(values, dtype=bool)
            for r in rows_c:
                values[bs.index(r["batch"]), col[r["marker"]]] = r["emd"]
                flagged[bs.index(r["batch"]), col[r["marker"]]] = r["flagged"]
            marks = {}
            for (b, m), info in missing.items():
                if m in col and b in bs:
                    flagged[bs.index(b), col[m]] = True
                    marks[(bs.index(b), col[m])] = GLYPH[info["kind"]]
            panels.append((f"{params['CONDITION_COLUMN']} = {c}", bs, markers, values, flagged, marks))
        paths.append(F.heatmap(panels, f"Each batch vs the other batches, within condition · {params['LABEL']}",
                               vmax, out / "heatmap_batches.png"))

    # Ridge plots for the most-flagged markers (then the largest distances).
    marker_flags = Counter(f["marker"] for f in summary["flags"] if f["level"] == "marker" and f["comparison"] != "missing")
    max_emd = defaultdict(float)
    for r in ref_rows + batch_rows:
        max_emd[r["marker"]] = max(max_emd[r["marker"]], r["emd"])
    ranked = sorted([m for m in markers if m in max_emd], key=lambda m: (-marker_flags.get(m, 0), -max_emd.get(m, 0.0)))
    for m in ranked[: max(0, params["MAX_PLOTS"])]:
        flagged_b = {f["batch"] for f in summary["flags"] if f["marker"] == m}
        series, pooled = {}, {}
        if ref_of:
            for b in batches:
                hs = [h for h in hists if h.batch == b and h.sample_id in ref_of]
                d = A.pooled_by(hs, m)
                if d is not None:
                    series[b] = {"reference sample": d}
            p = A.pool(s["reference sample"] for s in series.values())
            if p is not None:
                pooled["all batches"] = p
            what = "reference samples"
        else:
            conditions = sorted({h.condition for h in hists})
            for b in batches:
                per = {}
                for c in conditions:
                    d = A.pooled_by([h for h in hists if h.batch == b and h.condition == c], m)
                    if d is not None:
                        per[c] = d
                if per:
                    series[b] = per
            if len(conditions) == 1:
                p = A.pool(s[conditions[0]] for s in series.values() if conditions[0] in s)
                if p is not None:
                    pooled["all batches"] = p
            what = f"all cells by {params['CONDITION_COLUMN']}"
        if not series:
            continue
        name = f"density_{F.safe_name(m)}.png"
        paths.append(F.ridge(m, batches, series, pooled, flagged_b,
                             f"{m}: {what} by batch · {params['LABEL']}", out / name))
    return paths


if __name__ == "__main__":
    sys.exit(main())
