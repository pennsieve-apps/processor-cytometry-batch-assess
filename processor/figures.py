"""
PNG figures: batch x marker EMD heatmaps and per-marker ridge plots of the
distributions by batch. Colours follow one system: a single-hue sequential blue
ramp for EMD magnitude, categorical slots (fixed order) for conditions, the
reserved 'critical' red only for flags (always with a text mark, never colour
alone), muted inks for chrome.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from processor.assess import BIN_WIDTH, CENTERS

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
CRITICAL = "#d03b3b"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQUENTIAL = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": INK, "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold",
    })
    return plt


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", text).strip("_") or "marker"


def heatmap(panels: list[tuple[str, list[str], list[str], np.ndarray, np.ndarray, dict]], title: str, vmax: float,
            path: Path) -> Path:
    """panels: (subtitle, row labels, marker labels, emd matrix with NaN = no value, flagged bool matrix,
    {(row, col): glyph} for flagged cells without a distance, e.g. a missing marker)."""
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("emd", SEQUENTIAL)
    cmap.set_bad("#efeee9")
    n_cols = max(len(p[2]) for p in panels)
    heights = [max(len(p[1]), 1) for p in panels]
    fig_h = 1.3 + 0.32 * sum(heights) + 0.55 * len(panels)
    fig, axes = plt.subplots(len(panels), 1, figsize=(max(6.5, 0.36 * n_cols + 2.6), fig_h),
                             gridspec_kw={"height_ratios": heights}, squeeze=False)
    image = None
    for ax, (subtitle, rows, cols, values, flagged, marks) in zip(axes[:, 0], panels):
        image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=0, vmax=vmax, aspect="auto",
                          interpolation="nearest")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8, color=INK_2)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows, fontsize=8, color=INK_2)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(rows)), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=2)
        ax.tick_params(which="minor", length=0)
        if subtitle:
            ax.set_title(subtitle, loc="left", fontsize=9, fontweight="normal", color=INK_2)
        for i, j in zip(*np.nonzero(flagged)):
            ax.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, fill=False, edgecolor=CRITICAL, linewidth=1.6))
            v = values[i, j]
            ax.text(j, i, marks.get((i, j), "") if np.isnan(v) else f"{v:.2f}", ha="center", va="center", fontsize=6.5,
                    color=INK if np.isnan(v) or v < 0.55 * vmax else SURFACE, fontweight="bold")
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.02, 0.9, 0.97))
    cax = fig.add_axes((0.915, 0.25, 0.015, 0.5))
    cb = fig.colorbar(image, cax=cax)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=MUTED, length=0)
    cb.set_label("EMD (transformed units)", fontsize=8, color=INK_2)
    fig.text(0.01, 0.005, "Red outline: flagged (above median + K·MAD of the batches and MIN_EMD); ✕ = marker missing in this batch; ◆ = marker only in these batches (panel differs); grey = no value.",
             fontsize=7, color=MUTED)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _smooth(dist: np.ndarray, sigma_bins: float = 8.0) -> np.ndarray:
    half = int(4 * sigma_bins)
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / sigma_bins) ** 2)
    k /= k.sum()
    return np.convolve(dist, k, mode="same") / BIN_WIDTH


def ridge(marker: str, batches: list[str], series: dict[str, dict[str, np.ndarray]], pooled: dict[str, np.ndarray],
          flagged_batches: set[str], title: str, path: Path) -> Path:
    """series[batch][line name] = distribution; pooled[line name] = the all-batch reference drawn dashed in every row."""
    plt = _plt()
    line_names = sorted({n for b in series.values() for n in b})
    colours = {n: SERIES[i % len(SERIES)] for i, n in enumerate(line_names)}
    every = [d for b in series.values() for d in b.values()] + list(pooled.values())
    total = np.sum(every, axis=0)
    cdf = np.cumsum(total) / total.sum()
    lo = CENTERS[np.searchsorted(cdf, 0.002)]
    hi = CENTERS[min(np.searchsorted(cdf, 0.998), len(CENTERS) - 1)]
    pad = 0.08 * (hi - lo or 1.0)
    window = (CENTERS >= lo - pad) & (CENTERS <= hi + pad)
    x = CENTERS[window]
    curves = {b: {n: _smooth(d)[window] for n, d in s.items()} for b, s in series.items()}
    ref_curves = {n: _smooth(d)[window] for n, d in pooled.items()}
    peak = max([c.max() for s in curves.values() for c in s.values()] + [c.max() for c in ref_curves.values()] + [1e-9])
    step = 0.75
    fig, ax = plt.subplots(figsize=(6.5, 1.3 + 0.5 * len(batches)))
    for row, b in enumerate(reversed(batches)):
        base = row * step
        ax.axhline(base, color=GRID, linewidth=0.8, zorder=1)
        for n, c in curves.get(b, {}).items():
            y = base + c / peak * step * 1.6
            ax.fill_between(x, base, y, color=colours[n], alpha=0.18, linewidth=0, zorder=2 + row)
            ax.plot(x, y, color=colours[n], linewidth=1.6, zorder=3 + row)
        for n, c in ref_curves.items():
            ax.plot(x, base + c / peak * step * 1.6, color=INK_2, linewidth=0.9, linestyle=(0, (3, 2)), zorder=4 + row)
        is_flagged = b in flagged_batches
        ax.text(x[0] - 0.02 * (x[-1] - x[0]), base + 0.15 * step, ("▲ " if is_flagged else "") + b,
                ha="right", va="bottom", fontsize=8, color=CRITICAL if is_flagged else INK_2,
                fontweight="bold" if is_flagged else "normal")
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-0.1, (len(batches) - 1) * step + step * 1.9)
    ax.set_xlabel(f"{marker} (transformed)", color=INK_2)
    handles = [plt.Line2D([], [], color=colours[n], linewidth=1.6, label=n) for n in line_names]
    handles.append(plt.Line2D([], [], color=INK_2, linewidth=0.9, linestyle=(0, (3, 2)), label="all batches pooled"))
    handles.append(plt.Line2D([], [], color=CRITICAL, marker="^", linestyle="none", markersize=5, label="flagged batch"))
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=7, frameon=False,
              ncol=min(len(handles), 4), borderaxespad=0.2)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.subplots_adjust(left=0.2)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
