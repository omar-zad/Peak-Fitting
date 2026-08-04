#!/usr/bin/env python3
"""Plot a small set of representative GIWAXS frames without fitting.

The input profiles are the monitor-normalized FR/IP/OOP text files produced by
the preparation workflow.  The plots deliberately use the unsmoothed values:
one figure preserves the true intensity scale, while the second applies only a
constant vertical display offset so weak and overlapping features are visible.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault(
    "MPLCONFIGDIR", str(PROJECT_ROOT / ".cache" / "matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CUTS = ("FR", "IP", "OOP")


def parse_frames(value: str) -> list[int]:
    frames = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not frames:
        raise argparse.ArgumentTypeError("Provide at least one frame number.")
    if len(frames) != len(set(frames)):
        raise argparse.ArgumentTypeError("Frame numbers must be unique.")
    return frames


def locate_profile(processed_dir: Path, scan_id: str, cut: str) -> Path:
    matches = sorted(
        processed_dir.glob(f"scan{scan_id}_{cut}_*frames_norm.txt")
    )
    if len(matches) != 1:
        found = ", ".join(path.name for path in matches) or "none"
        raise FileNotFoundError(
            f"Expected one normalized {cut} profile for scan {scan_id}; "
            f"found {found}."
        )
    return matches[0]


def load_profile(path: Path, frames: list[int]) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    table = pd.read_csv(path, sep="\t")
    q = pd.to_numeric(table.iloc[:, 0], errors="coerce").to_numpy(float)
    profiles: dict[int, np.ndarray] = {}
    for frame in frames:
        column = str(frame)
        if column not in table.columns:
            raise KeyError(f"Frame {frame} is missing from {path.name}.")
        profiles[frame] = pd.to_numeric(
            table[column], errors="coerce"
        ).to_numpy(float)
    return q, profiles


def load_manifest(
    processed_dir: Path,
    scan_id: str,
    frames: list[int],
) -> tuple[Path, pd.DataFrame]:
    path = processed_dir / f"scan{scan_id}_frame_manifest.csv"
    manifest = pd.read_csv(path)
    manifest["frame"] = pd.to_numeric(manifest["frame"], errors="raise").astype(int)
    selected = manifest.set_index("frame").reindex(frames)
    if selected.isna().all(axis=1).any():
        missing = selected.index[selected.isna().all(axis=1)].tolist()
        raise KeyError(f"Frames missing from manifest: {missing}")
    return path, selected.reset_index()


def frame_label(row: pd.Series) -> str:
    frame = int(row["frame"])
    temperature = float(row["measured_temperature_C"])
    high_time = pd.to_numeric(
        pd.Series([row.get("time_at_high_temperature_s")]),
        errors="coerce",
    ).iloc[0]
    if np.isfinite(high_time):
        time_text = f"HT+{high_time / 60:.1f} min"
    else:
        acquisition = float(row["acquisition_time_s"])
        time_text = f"t={acquisition / 60:.1f} min"
    return f"f{frame} | {temperature:.1f} °C | {time_text}"


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 10,
            "legend.fontsize": 8.2,
            "figure.titlesize": 15,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
        }
    )


def finite_xy(
    q: np.ndarray,
    intensity: np.ndarray,
    q_limits: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    keep = (
        np.isfinite(q)
        & np.isfinite(intensity)
        & (q >= q_limits[0])
        & (q <= q_limits[1])
    )
    return q[keep], intensity[keep]


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_true_intensity(
    scan_id: str,
    sample_label: str,
    data: dict[str, tuple[np.ndarray, dict[int, np.ndarray]]],
    frames: list[int],
    labels: dict[int, str],
    colors: dict[int, tuple[float, float, float, float]],
    output_stem: Path,
) -> None:
    ranges = [
        ("Full q range", (0.0, 2.5)),
        ("Low-q detail", (0.22, 0.60)),
        ("High-q detail", (1.42, 2.20)),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(17.2, 11.2))

    for row, cut in enumerate(CUTS):
        q, profiles = data[cut]
        for column, (title, limits) in enumerate(ranges):
            ax = axes[row, column]
            for frame in frames:
                x, y = finite_xy(q, profiles[frame], limits)
                ax.plot(
                    x,
                    y,
                    color=colors[frame],
                    lw=1.25,
                    alpha=0.94,
                    label=labels[frame],
                )
            if row == 0:
                ax.set_title(title)
            if row == 2:
                ax.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
            if column == 0:
                ax.set_ylabel(f"{cut} intensity (a.u.)")
            ax.set_xlim(*limits)
            ax.grid(alpha=0.18, lw=0.55)

    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.945),
    )
    panel_letters = iter("abcdefghi")
    for ax in axes.flat:
        ax.text(
            0.015,
            0.965,
            f"({next(panel_letters)})",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontweight="bold",
        )

    fig.suptitle(
        f"{sample_label} (scan {scan_id}): monitor-normalized GIWAXS profiles",
        y=0.995,
    )
    fig.text(
        0.5,
        0.965,
        "Unsmoothed; d5i-monitor normalized; no vertical offsets",
        ha="center",
        va="top",
        color="#555555",
    )
    fig.subplots_adjust(
        left=0.07,
        right=0.985,
        bottom=0.07,
        top=0.865,
        hspace=0.30,
        wspace=0.20,
    )
    save_figure(fig, output_stem)


def robust_offset_step(
    q: np.ndarray,
    profiles: dict[int, np.ndarray],
    frames: list[int],
    q_limits: tuple[float, float],
) -> float:
    values: list[np.ndarray] = []
    for frame in frames:
        _, y = finite_xy(q, profiles[frame], q_limits)
        values.append(y)
    combined = np.concatenate(values)
    low, high = np.nanpercentile(combined, [5, 99])
    span = max(float(high - low), np.finfo(float).eps)
    return 0.34 * span


def plot_stacked(
    scan_id: str,
    sample_label: str,
    data: dict[str, tuple[np.ndarray, dict[int, np.ndarray]]],
    frames: list[int],
    labels: dict[int, str],
    colors: dict[int, tuple[float, float, float, float]],
    output_stem: Path,
) -> None:
    ranges = [
        ("Full q range", (0.0, 2.5)),
        ("Low-q detail", (0.22, 0.60)),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(15.5, 11.0))

    for row, cut in enumerate(CUTS):
        q, profiles = data[cut]
        full_step = robust_offset_step(q, profiles, frames, ranges[0][1])
        low_step = robust_offset_step(q, profiles, frames, ranges[1][1])
        for column, ((title, limits), step) in enumerate(
            zip(ranges, (full_step, low_step))
        ):
            ax = axes[row, column]
            for offset_index, frame in enumerate(frames):
                x, y = finite_xy(q, profiles[frame], limits)
                ax.plot(
                    x,
                    y + offset_index * step,
                    color=colors[frame],
                    lw=1.25,
                    alpha=0.96,
                    label=labels[frame],
                )
            if row == 0:
                ax.set_title(title)
            if row == 2:
                ax.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
            ax.set_ylabel(f"{cut} intensity + offset")
            ax.set_xlim(*limits)
            ax.grid(alpha=0.18, lw=0.55)

    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.945),
    )
    panel_letters = iter("abcdef")
    for ax in axes.flat:
        ax.text(
            0.015,
            0.965,
            f"({next(panel_letters)})",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontweight="bold",
        )

    fig.suptitle(
        f"{sample_label} (scan {scan_id}): GIWAXS profile evolution",
        y=0.995,
    )
    fig.text(
        0.5,
        0.965,
        "Unsmoothed and d5i-monitor normalized; offsets aid visibility only",
        ha="center",
        va="top",
        color="#555555",
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.07,
        top=0.865,
        hspace=0.30,
        wspace=0.18,
    )
    save_figure(fig, output_stem)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot representative unsmoothed FR/IP/OOP frames."
    )
    parser.add_argument(
        "--processed-dir",
        type=Path,
        required=True,
        help="Folder containing normalized profiles and frame manifest.",
    )
    parser.add_argument("--scan-id", required=True)
    parser.add_argument(
        "--sample-label",
        default=None,
        help="Sample name displayed in the dissertation-draft figures.",
    )
    parser.add_argument(
        "--frames",
        type=parse_frames,
        required=True,
        help="Comma-separated zero-based frames, e.g. 0,25,50,57,71,85,99.",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    sample_label = args.sample_label or f"Scan {args.scan_id}"

    processed_dir = args.processed_dir.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_path, selected_manifest = load_manifest(
        processed_dir,
        args.scan_id,
        args.frames,
    )
    labels = {
        int(row["frame"]): frame_label(row)
        for _, row in selected_manifest.iterrows()
    }
    palette = plt.get_cmap("viridis")(
        np.linspace(0.08, 0.92, len(args.frames))
    )
    colors = {
        frame: tuple(palette[index])
        for index, frame in enumerate(args.frames)
    }

    profile_paths = {
        cut: locate_profile(processed_dir, args.scan_id, cut)
        for cut in CUTS
    }
    data = {
        cut: load_profile(path, args.frames)
        for cut, path in profile_paths.items()
    }

    apply_style()
    plot_true_intensity(
        args.scan_id,
        sample_label,
        data,
        args.frames,
        labels,
        colors,
        outdir / "key_frames_true_intensity",
    )
    plot_stacked(
        args.scan_id,
        sample_label,
        data,
        args.frames,
        labels,
        colors,
        outdir / "key_frames_stacked",
    )

    selected_manifest.to_csv(
        outdir / "selected_frame_manifest.csv",
        index=False,
    )
    provenance = {
        "scan_id": args.scan_id,
        "sample_label": sample_label,
        "frames": args.frames,
        "profiles": {
            cut: str(path) for cut, path in profile_paths.items()
        },
        "frame_manifest": str(manifest_path),
        "processing": (
            "Unsmoothed d5i-monitor-normalized profiles. "
            "The true-intensity figure has no offsets. The stacked figure "
            "adds constant display offsets only."
        ),
    }
    (outdir / "plot_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n",
        encoding="utf-8",
    )

    frame_groups: list[tuple[str, list[int]]] = []
    if "phase" in selected_manifest.columns:
        for phase in selected_manifest["phase"].dropna().astype(str).unique():
            group = selected_manifest.loc[
                selected_manifest["phase"].astype(str) == phase,
                "frame",
            ].astype(int).tolist()
            if group:
                frame_groups.append((phase, group))
    if not frame_groups:
        frame_groups = [("the measured series", args.frames)]
    phase_sentence = "; ".join(
        f"frames {frames} represent {phase}"
        for phase, frames in frame_groups
    )
    captions = (
        "# Dissertation-draft figure captions\n\n"
        "## Quantitative, un-offset profiles\n\n"
        f"Representative monitor-normalized one-dimensional GIWAXS profiles "
        f"for {sample_label} (scan {args.scan_id}) in the full-radial (FR), "
        "in-plane (IP), and out-of-plane (OOP) integration sectors. "
        f"The temperature-trace phase labels are: {phase_sentence}. "
        "These labels should be checked against the experimental protocol "
        "before final dissertation use. The profiles are unsmoothed and "
        "normalized using the incident-beam d5i monitor. No vertical offsets "
        "have been applied. Full-q profiles are accompanied by expanded "
        "low- and high-q regions to show changes in peak presence and "
        "position.\n\n"
        "## Vertically offset profiles\n\n"
        f"Evolution of representative monitor-normalized GIWAXS profiles for "
        f"{sample_label} (scan {args.scan_id}) in the FR, IP, and OOP "
        "integration sectors. Constant vertical offsets have been applied "
        "for visual clarity only; the profiles were not smoothed or "
        "individually rescaled. The offset presentation highlights the "
        "appearance, disappearance, splitting, and displacement of weak "
        "features across the selected experimental phases.\n\n"
        "Draft note: replace the automatic phase labels with the final agreed "
        "experimental wording wherever the beamline log or protocol provides "
        "a more specific description.\n"
    )
    (outdir / "dissertation_figure_captions.md").write_text(
        captions,
        encoding="utf-8",
    )

    print(f"Saved key-frame plots to {outdir}")


if __name__ == "__main__":
    main()
