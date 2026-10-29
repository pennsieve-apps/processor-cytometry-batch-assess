"""
The batch-effect measures.

Each sample is read once, subsampled (seeded) and reduced to one histogram per
marker on a fixed fine grid over the transformed scale. Distributions are then
pooled (every sample, and every batch, weighs the same) and compared with the
1-D Earth Mover's Distance (scipy.stats.wasserstein_distance on the bin centres,
weighted by the bin masses). Binning moves each value by at most half a bin, so
the EMD is exact to within one bin width (BIN_WIDTH = 0.005 transformed units).

Two comparisons:

- reference: when reference samples are named (REFERENCE_COLUMN /
  REFERENCE_VALUES, e.g. a replicate of the same donor aliquot in every batch),
  each batch's reference sample(s) against the reference samples of all
  batches pooled. Same cells in every batch, so any distance is technical.
- batches: always; within each condition, each batch's cells (its samples
  pooled) against the cells of all *other* batches pooled, so a difference
  between conditions (biology) isn't counted as batch effect.

A batch x marker distance is flagged when it is above a robust threshold over
the batches of the same comparison (median + K * MAD, MAD scaled by 1.4826 so it
estimates the SD for normal data) and above MIN_EMD.
"""
from __future__ import annotations

import math
import zlib
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
from scipy.stats import wasserstein_distance

GRID_LO, GRID_HI, BIN_WIDTH = -10.0, 20.0, 0.005
N_BINS = int(round((GRID_HI - GRID_LO) / BIN_WIDTH))
CENTERS = GRID_LO + BIN_WIDTH * (np.arange(N_BINS) + 0.5)
MAD_SCALE = 1.4826


@dataclass
class SampleHist:
    sample_id: str
    batch: str
    condition: str
    counts: dict[str, np.ndarray]        # marker -> int64 counts over the grid
    n_events: int
    n_cells: int                         # events used (after subsampling)
    clipped: dict[str, float] = field(default_factory=dict)   # marker -> fraction outside the grid
    absent: list[str] = field(default_factory=list)           # panel markers this file lacks


def sample_seed(seed: int, sample_id: str) -> int:
    """A per-sample seed, so a sample's subsample doesn't depend on which other samples are in the run."""
    return (int(seed) * 1_000_003 + zlib.crc32(sample_id.encode())) % (2**32)


def subsample_rows(n: int, cells: int, rng: np.random.Generator) -> np.ndarray | None:
    if cells <= 0 or n <= cells:
        return None
    return np.sort(rng.choice(n, size=cells, replace=False))


def histogram(values: np.ndarray) -> tuple[np.ndarray, float]:
    finite = values[np.isfinite(values)]
    idx = np.floor((finite.astype(np.float64) - GRID_LO) / BIN_WIDTH).astype(np.int64)
    outside = np.count_nonzero((idx < 0) | (idx >= N_BINS))
    np.clip(idx, 0, N_BINS - 1, out=idx)
    return np.bincount(idx, minlength=N_BINS), (outside / finite.size if finite.size else 0.0)


def sample_histograms(sample_id: str, batch: str, condition: str, pnn: list[str], events: np.ndarray,
                      markers: list[tuple[str, str]], cells: int, seed: int) -> SampleHist:
    rng = np.random.default_rng(sample_seed(seed, sample_id))
    rows = subsample_rows(events.shape[0], cells, rng)
    used = events if rows is None else events[rows]
    col = {name: i for i, name in enumerate(pnn)}
    counts, clipped, absent = {}, {}, []
    for channel, marker in markers:
        if channel not in col:
            absent.append(marker)
            continue
        counts[marker], clipped[marker] = histogram(used[:, col[channel]])
    return SampleHist(sample_id, batch, condition, counts, events.shape[0], used.shape[0], clipped, absent)


def normalised(counts: np.ndarray) -> np.ndarray | None:
    total = counts.sum()
    return counts / total if total > 0 else None


def pool(dists: Iterable[np.ndarray | None]) -> np.ndarray | None:
    """Equal-weight mixture of the given distributions (each sums to 1)."""
    dists = [d for d in dists if d is not None]
    if not dists:
        return None
    return np.mean(dists, axis=0)


def emd(p: np.ndarray, q: np.ndarray) -> float:
    support = (p > 0) | (q > 0)
    return float(wasserstein_distance(CENTERS[support], CENTERS[support], p[support], q[support]))


def pooled_by(hists: list[SampleHist], marker: str) -> np.ndarray | None:
    return pool(normalised(h.counts[marker]) for h in hists if marker in h.counts)


def compare_reference(hists: list[SampleHist], ref_of: dict[str, str], markers: list[str]) -> list[dict]:
    """Each batch's reference sample(s) vs the reference samples of all batches, per reference value."""
    groups: dict[str, dict[str, list[SampleHist]]] = defaultdict(lambda: defaultdict(list))
    for h in hists:
        if h.sample_id in ref_of:
            groups[ref_of[h.sample_id]][h.batch].append(h)
    rows = []
    for ref, by_batch in sorted(groups.items()):
        for marker in markers:
            per_batch = {b: pooled_by(hs, marker) for b, hs in by_batch.items()}
            per_batch = {b: d for b, d in per_batch.items() if d is not None}
            if len(per_batch) < 2:
                continue
            reference = pool(per_batch.values())
            for b in sorted(per_batch):
                hs = [h for h in by_batch[b] if marker in h.counts]
                rows.append({"batch": b, "reference": ref, "marker": marker, "emd": emd(per_batch[b], reference),
                             "n_samples": len(hs), "n_cells": int(sum(h.n_cells for h in hs))})
    return rows


def compare_batches(hists: list[SampleHist], markers: list[str]) -> list[dict]:
    """Within each condition: each batch's pooled cells vs all other batches' pooled cells."""
    groups: dict[str, dict[str, list[SampleHist]]] = defaultdict(lambda: defaultdict(list))
    for h in hists:
        groups[h.condition][h.batch].append(h)
    rows = []
    for condition, by_batch in sorted(groups.items()):
        for marker in markers:
            per_batch = {b: pooled_by(hs, marker) for b, hs in by_batch.items()}
            per_batch = {b: d for b, d in per_batch.items() if d is not None}
            if len(per_batch) < 2:
                continue
            for b in sorted(per_batch):
                others = pool(d for o, d in per_batch.items() if o != b)
                hs = [h for h in by_batch[b] if marker in h.counts]
                rows.append({"batch": b, "condition": condition, "marker": marker, "emd": emd(per_batch[b], others),
                             "n_samples": len(hs), "n_cells": int(sum(h.n_cells for h in hs))})
    return rows


def robust_threshold(values: list[float], k: float) -> tuple[float, float, float]:
    """(median, scaled MAD, median + k * MAD)."""
    arr = np.asarray(values, dtype=float)
    med = float(np.median(arr))
    mad = float(MAD_SCALE * np.median(np.abs(arr - med)))
    return med, mad, med + k * mad


def mark_outliers(rows: list[dict], group_keys: tuple[str, ...], k: float, min_emd: float, min_batches: int = 3) -> None:
    """Adds threshold / flagged to each row; the threshold is over the batches of the same group."""
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        groups[tuple(r[g] for g in group_keys)].append(r)
    for members in groups.values():
        if len(members) < min_batches:
            for r in members:
                r["threshold"], r["flagged"] = float("nan"), False
            continue
        _, _, thr = robust_threshold([r["emd"] for r in members], k)
        for r in members:
            r["threshold"] = thr
            r["flagged"] = bool(r["emd"] > thr and r["emd"] > min_emd)


def missing_by_batch(hists: list[SampleHist], prep_missing: dict[str, set[str]],
                     batch_of: dict[str, str]) -> dict[tuple[str, str], dict]:
    """(batch, marker) -> {kind, n_missing, n_samples, source, reason} for markers some batches lack.

    A marker lacking in at most half of the batches flags those batches (missing); one lacking in
    most batches flags the batches that have it (only_here). Lacking in all batches is not a flag.
    """
    n_in_batch: dict[str, int] = defaultdict(int)
    for sid, b in batch_of.items():
        n_in_batch[b] += 1
    lacking: dict[tuple[str, str], set[str]] = defaultdict(set)
    sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    for h in hists:
        for m in h.absent:
            lacking[(h.batch, m)].add(h.sample_id)
            sources[(h.batch, m)].add("fcs")
    for marker, sids in prep_missing.items():
        for sid in sids:
            if sid in batch_of:
                lacking[(batch_of[sid], marker)].add(sid)
                sources[(batch_of[sid], marker)].add("prep.json")
    batches = set(batch_of.values())
    by_marker: dict[str, set[str]] = defaultdict(set)
    for (b, m) in lacking:
        if len(lacking[(b, m)]) and b in batches:
            by_marker[m].add(b)
    out = {}
    for m, lacking_batches in by_marker.items():
        if lacking_batches >= batches:
            continue   # missing from every batch: a panel difference, not a batch effect
        having = batches - lacking_batches
        if len(lacking_batches) <= len(having):
            for b in lacking_batches:
                sids = lacking[(b, m)]
                out[(b, m)] = {"kind": "missing", "n_missing": len(sids), "n_samples": n_in_batch[b],
                               "source": "+".join(sorted(sources[(b, m)])),
                               "reason": f"marker missing in {len(sids)} of {n_in_batch[b]} samples "
                                         f"({'+'.join(sorted(sources[(b, m)]))})"}
        else:
            # Most batches lack it: the batches that have it are the odd ones (a panel change).
            for b in having:
                out[(b, m)] = {"kind": "only_here", "n_missing": 0, "n_samples": n_in_batch[b],
                               "source": "+".join(sorted(set().union(*(sources[(o, m)] for o in lacking_batches)))),
                               "reason": f"marker only in {len(having)} of {len(batches)} batches (panel differs)"}
    return out


def summarise(label: str, batches: list[str], markers: list[str], ref_rows: list[dict], batch_rows: list[dict],
              missing: dict[tuple[str, str], dict], batch_fraction: float) -> tuple[list[dict], dict]:
    """flags.csv rows and the summary.json body."""
    flags: list[dict] = []
    for r in ref_rows:
        if r["flagged"]:
            flags.append({"level": "marker", "comparison": "reference", "batch": r["batch"], "condition": "",
                          "reference": r["reference"], "marker": r["marker"], "emd": r["emd"], "threshold": r["threshold"],
                          "reason": "reference sample differs from the pooled reference samples"})
    for r in batch_rows:
        if r["flagged"]:
            flags.append({"level": "marker", "comparison": "batches", "batch": r["batch"], "condition": r["condition"],
                          "reference": "", "marker": r["marker"], "emd": r["emd"], "threshold": r["threshold"],
                          "reason": "batch differs from the other batches within the condition"})
    for (b, m), info in sorted(missing.items()):
        flags.append({"level": "marker", "comparison": "missing", "batch": b, "condition": "", "reference": "",
                      "marker": m, "emd": float("nan"), "threshold": float("nan"), "reason": info["reason"]})

    def mean_by_batch(rows):
        acc: dict[str, list[float]] = defaultdict(list)
        for r in rows:
            acc[r["batch"]].append(r["emd"])
        return {b: float(np.mean(v)) for b, v in acc.items()}

    ref_mean, batch_mean = mean_by_batch(ref_rows), mean_by_batch(batch_rows)
    n_markers = len(markers)
    need = max(1, math.ceil(batch_fraction * n_markers)) if n_markers else 1
    per_batch = {}
    for b in batches:
        flagged = sorted({f["marker"] for f in flags if f["batch"] == b})
        score = ref_mean.get(b, batch_mean.get(b))
        entry = {
            "score": score,
            "score_source": "reference" if b in ref_mean else ("batches" if b in batch_mean else None),
            "mean_emd_reference": ref_mean.get(b),
            "mean_emd_batches": batch_mean.get(b),
            "n_flagged_markers": len(flagged),
            "flagged_markers": flagged,
            "flagged": len(flagged) >= need,
        }
        per_batch[b] = entry
        if entry["flagged"]:
            flags.append({"level": "batch", "comparison": "summary", "batch": b, "condition": "", "reference": "",
                          "marker": "", "emd": score if score is not None else float("nan"), "threshold": float("nan"),
                          "reason": f"{len(flagged)} of {n_markers} markers flagged (>= {need})"})

    def stats(rows):
        if not rows:
            return None
        v = np.array([r["emd"] for r in rows])
        return {"n": int(v.size), "median_emd": float(np.median(v)), "mean_emd": float(v.mean()),
                "max_emd": float(v.max()), "n_flagged": int(sum(r["flagged"] for r in rows))}

    ref_stats, batch_stats = stats(ref_rows), stats(batch_rows)
    headline = ref_stats or batch_stats
    summary = {
        "label": label,
        "n_batches": len(batches),
        "n_markers": n_markers,
        "median_emd": headline["median_emd"] if headline else None,
        "median_emd_source": "reference" if ref_stats else ("batches" if batch_stats else None),
        "comparisons": {"reference": ref_stats, "batches": batch_stats},
        "per_batch": per_batch,
        "flagged_batches": [b for b, e in per_batch.items() if e["flagged"]],
        "flags": flags,
    }
    return flags, summary
