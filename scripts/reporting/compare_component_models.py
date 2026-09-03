#!/usr/bin/env python3
"""Compare saved component-count runs for one common fitting window."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".cache" / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_run(specification: str) -> tuple[str, Path]:
    if "=" not in specification:
        raise argparse.ArgumentTypeError("Expected LABEL=RESULT_DIRECTORY")
    label, path_text = specification.split("=", 1)
    return label.strip(), Path(path_text).expanduser()


def selected_row(table: pd.DataFrame, window: str) -> pd.Series:
    rows = table.loc[table["window"].astype(str).eq(window)].copy()
    if rows.empty:
        raise ValueError(f"No model-selection row found for window {window!r}")
    selected = rows.loc[
        rows["selected"].astype(str).str.lower().isin(("true", "1", "yes"))
    ]
    return (selected if not selected.empty else rows).sort_values("bic").iloc[0]


def collect(label: str, directory: Path, window: str) -> dict[str, object]:
    models = pd.read_csv(directory / "model_selection.csv")
    peaks = pd.read_csv(directory / "peak_parameters.csv")
    curves = pd.read_csv(directory / "fitted_curves.csv")
    model = selected_row(models, window)
    peaks = peaks.loc[peaks["window"].astype(str).eq(window)]
    curves = curves.loc[curves["window"].astype(str).eq(window)]
    return {
        "label": label,
        "result_directory": str(directory),
        "n_components": int(peaks["peak"].nunique()),
        "n_frames": int(peaks["frame"].nunique()),
        "profile_family": str(model["profile"]),
        "bic": float(model["bic"]),
        "n_points": int(model["n_points"]),
        "n_parameters": int(model["n_parameters"]),
        "q_min_Ainv": float(curves["q_Ainv"].min()),
        "q_max_Ainv": float(curves["q_Ainv"].max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True, type=parse_run)
    parser.add_argument("--window", default="lamellar_low_q")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--title", default="Component-count model comparison")
    args = parser.parse_args()

    records = [collect(label, path, args.window) for label, path in args.run]
    table = pd.DataFrame.from_records(records).sort_values("n_components")
    signatures = table[["n_points", "q_min_Ainv", "q_max_Ainv"]].round(8)
    if len(signatures.drop_duplicates()) != 1:
        raise ValueError(
            "Runs are not comparable: fitted point counts or q-ranges differ"
        )

    table["delta_bic"] = table["bic"] - table["bic"].min()
    table["preferred"] = table["delta_bic"].eq(0)
    args.outdir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.outdir / "component_model_comparison.csv", index=False)

    figure, axis = plt.subplots(figsize=(7.2, 4.7), constrained_layout=True)
    colors = np.where(table["preferred"], "#176B87", "#A9B8C6")
    bars = axis.bar(table["n_components"], table["delta_bic"], color=colors)
    for bar, row in zip(bars, table.itertuples()):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"Delta BIC={row.delta_bic:.1f}\nBIC={row.bic:.1f}\n{row.profile_family}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    axis.axhline(6, color="#B3413E", linestyle="--", linewidth=1, label="Delta BIC = 6")
    axis.set_xticks(table["n_components"])
    axis.set_xlabel("Number of fitted peak components")
    axis.set_ylabel("Delta BIC (0 is preferred)")
    maximum_delta = float(table["delta_bic"].max())
    axis.set_ylim(0, maximum_delta * 1.20 if maximum_delta > 0 else 1)
    preferred = table.loc[table["preferred"]].iloc[0]
    axis.scatter(
        [preferred["n_components"]],
        [0],
        color="#176B87",
        marker="D",
        s=34,
        zorder=4,
        clip_on=False,
    )
    axis.set_title(args.title, pad=12)
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    figure.savefig(args.outdir / "component_model_comparison.png", dpi=220)
    figure.savefig(args.outdir / "component_model_comparison.pdf")
    plt.close(figure)


if __name__ == "__main__":
    main()
