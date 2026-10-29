"""
Finding the inputs and reading the metadata.

Linked layout (placement: subfolder): INPUT_DIR/fcs/<sample_id>.fcs and
INPUT_DIR/metadata/{samples.csv, panel.csv, prep.json}. Without links the same
files are found by name anywhere under INPUT_DIR (symlinks followed).
"""
from __future__ import annotations

import json
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

log = logging.getLogger("processor-cytometry-batch-assess")


@dataclass
class Sample:
    sample_id: str
    path: Path
    batch: str
    condition: str
    row: dict[str, Any] = field(default_factory=dict)


@dataclass
class Inputs:
    samples: list[Sample]
    markers: list[tuple[str, str]]               # (channel $PnN, marker name), panel order
    missing: dict[str, set[str]]                 # marker -> sample_ids that lack it (prep.json)
    files_seen: list[Path]
    samples_csv: Path | None
    panel_csv: Path | None
    prep_json: Path | None
    skipped: dict[str, str]                      # sample_id or file -> reason


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinks (INPUT_DIR entries are symlinks), sorted."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return found


def _find(root: Path, port: str, predicate) -> list[Path]:
    """Files of a port: its own folder when it exists, else anywhere under INPUT_DIR."""
    folder = root / port
    base = folder if folder.is_dir() else root
    return [p for p in walk(base) if predicate(p)]


def _one(root: Path, name: str) -> Path | None:
    hits = _find(root, "metadata", lambda p: p.name.lower() == name)
    if len(hits) > 1:
        log.warning("%d files named %s; using %s", len(hits), name, hits[0])
    return hits[0] if hits else None


def _text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "NA"
    text = str(value).strip()
    return text or "NA"


def read_panel(path: Path) -> list[tuple[str, str]]:
    panel = pd.read_csv(path, dtype=str, keep_default_na=False)
    cols = {c.lower().strip(): c for c in panel.columns}
    if "channel" not in cols:
        raise ValueError(f"{path.name}: no 'channel' column")
    rows = []
    for _, r in panel.iterrows():
        kind = r[cols["type"]].strip().lower() if "type" in cols else "marker"
        if kind != "marker":
            continue
        channel = r[cols["channel"]].strip()
        marker = r[cols["marker"]].strip() if "marker" in cols else ""
        rows.append((channel, marker or channel))
    return rows


def read_prep_missing(path: Path, channel_to_marker: dict[str, str]) -> dict[str, set[str]]:
    """Markers that some samples lack, from fcs-prep's prep.json, as marker -> {sample_id lacking it}.

    fcs-prep keeps the channels every sample has with the same marker, and lists per sample the
    channels it *had* that were dropped: samples.<id>.channels_dropped = [{channel, marker, reason}].
    A dropped marker is missing from the samples that don't list it (a marker that moved to another
    channel, as CD28 does in the demo data, is listed by every sample and so is not missing).
    Entries named *missing* ([channel or {channel, marker}]) are taken as lacking, for other writers.
    """
    with open(path) as fh:
        data = json.load(fh)
    samples = data.get("samples") if isinstance(data, dict) else None
    if isinstance(samples, list):
        samples = {str(e.get("sample_id")): e for e in samples if isinstance(e, dict) and e.get("sample_id")}
    if not isinstance(samples, dict):
        return {}

    def name(item) -> str:
        if isinstance(item, dict):
            channel = str(item.get("channel") or "")
            return str(item.get("marker") or channel_to_marker.get(channel) or channel)
        return channel_to_marker.get(str(item), str(item))

    had: dict[str, set[str]] = {}
    lacking: dict[str, set[str]] = {}
    for sid, entry in samples.items():
        if not isinstance(entry, dict):
            continue
        for key, value in entry.items():
            if not isinstance(value, list):
                continue
            k = key.lower()
            if "dropped" in k:
                for item in value:
                    if name(item):
                        had.setdefault(name(item), set()).add(str(sid))
            elif "missing" in k:
                for item in value:
                    if name(item):
                        lacking.setdefault(name(item), set()).add(str(sid))
    every = {str(s) for s in samples}
    for marker, ids in had.items():
        if every - ids:
            lacking.setdefault(marker, set()).update(every - ids)
    return lacking


def discover(input_dir: Path, batch_column: str, condition_column: str) -> Inputs:
    fcs_files = _find(input_dir, "fcs", lambda p: p.suffix.lower() == ".fcs")
    samples_csv = _one(input_dir, "samples.csv")
    panel_csv = _one(input_dir, "panel.csv")
    prep_json = _one(input_dir, "prep.json")
    skipped: dict[str, str] = {}

    by_stem: dict[str, Path] = {}
    for p in fcs_files:
        if p.stem in by_stem:
            log.warning("two FCS files for sample %s; using %s", p.stem, by_stem[p.stem])
            continue
        by_stem[p.stem] = p

    samples: list[Sample] = []
    if samples_csv is not None:
        table = pd.read_csv(samples_csv, dtype=str, keep_default_na=False)
        if "sample_id" not in table.columns:
            raise ValueError(f"{samples_csv.name}: no sample_id column")
        if batch_column not in table.columns:
            raise ValueError(f"{samples_csv.name}: no {batch_column!r} column")
        listed = set()
        for _, r in table.iterrows():
            row = {k: _text(v) for k, v in r.items()}
            sid = row["sample_id"]
            listed.add(sid)
            path = by_stem.get(sid)
            if path is None and row.get("file") not in (None, "NA"):
                path = next((p for p in fcs_files if p.name == row["file"]), None)
            if path is None:
                skipped[sid] = "listed in samples.csv but no FCS file found"
                continue
            condition = row.get(condition_column, "NA") if condition_column else "NA"
            samples.append(Sample(sid, path, row[batch_column], condition, row))
        for stem, path in by_stem.items():
            if stem not in listed and not any(s.path == path for s in samples):
                skipped[stem] = "FCS file not listed in samples.csv"
    else:
        # No sample sheet: batch = the file's $DATE keyword, as fcs-prep does.
        from processor.fcs import read_header

        log.warning("no samples.csv: batch = $DATE keyword of each file, condition = NA")
        for stem, path in sorted(by_stem.items()):
            date = read_header(path).keywords.get("date") or "NA"
            samples.append(Sample(stem, path, date, "NA", {"sample_id": stem, "batch": date}))

    markers: list[tuple[str, str]] = []
    if panel_csv is not None:
        markers = read_panel(panel_csv)
    elif samples:
        # No panel: every channel with a $PnS label that isn't scatter or time.
        from processor.fcs import read_header

        log.warning("no panel.csv: using channels with a $PnS label (not scatter/time) of %s", samples[0].path.name)
        h = read_header(samples[0].path)
        for n, s in zip(h.pnn, h.pns):
            if s and not re.match(r"^(FSC|SSC|Time)", n, re.I):
                markers.append((n, s))

    missing: dict[str, set[str]] = {}
    if prep_json is not None:
        try:
            missing = read_prep_missing(prep_json, dict(markers))
        except (OSError, ValueError) as exc:
            log.warning("cannot read %s: %s", prep_json, exc)
    return Inputs(samples, markers, missing, fcs_files, samples_csv, panel_csv, prep_json, skipped)
