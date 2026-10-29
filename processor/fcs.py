"""
Reading FCS files with FlowIO: the channel names ($PnN, $PnS), a few keywords,
and the events as a float32 matrix (one row per event). Values are taken as
stored; fcs-prep has already arcsinh-transformed the marker channels and
downstream steps never re-transform (see the workflow conventions).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class FcsChannels:
    pnn: list[str]
    pns: list[str]
    keywords: dict[str, str]
    n_events: int


def _channels(flow) -> tuple[list[str], list[str]]:
    order = sorted(flow.channels, key=int)
    pnn = [str(flow.channels[k].get("PnN") or "") for k in order]
    pns = [str(flow.channels[k].get("PnS") or "") for k in order]
    return pnn, pns


def read_header(path: Path) -> FcsChannels:
    import flowio

    flow = flowio.FlowData(str(path), only_text=True)
    pnn, pns = _channels(flow)
    return FcsChannels(pnn, pns, dict(flow.text), int(flow.event_count))


def read_events(path: Path) -> tuple[FcsChannels, np.ndarray]:
    """Channels and the event matrix (n_events x n_channels, float32)."""
    import flowio

    flow = flowio.FlowData(str(path))
    pnn, pns = _channels(flow)
    events = np.asarray(flow.events, dtype=np.float32).reshape(-1, flow.channel_count)
    del flow
    return FcsChannels(pnn, pns, {}, events.shape[0]), events
