#!/usr/bin/env python3
"""Build a read-only, supervisor-ready QC package from GIWAXS fit outputs.

The script does not refit or alter any source result. It combines the final
per-cut CSV files, applies transparent reporting gates, compares available
low-q candidate models, and creates diagnostic figures plus a Markdown report.

Default final inputs:
    results/drop40/analysis_9frame_fast/FR_revised
    results/drop40/analysis_9frame_fast/IP_revised
    results/drop40/analysis_9frame_fast/OOP_revised

Example:
    python scripts/reporting/build_nine_frame_qc_report.py

Alternative inputs can be supplied repeatedly:
    python scripts/reporting/build_nine_frame_qc_report.py \
        --input FR=results/drop40/analysis_9frame_fast/FR_revised \
        --input IP=results/drop40/analysis_9frame_fast/IP_revised \
        --input OOP=results/drop40/analysis_9frame_fast/OOP_revised

Terminology
-----------
The optional value 2*pi/FWHM is labelled throughout as an *uncorrected
reciprocal-space apparent coherence length*. It is not labelled CCL because
instrumental broadening correction, peak assignment, and a manual
background/window sensitivity review remain necessary.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(PROJECT_ROOT / ".cache" / "matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TABLE_NAMES = (
    "model_selection",
    "fit_quality",
    "peak_parameters",
    "fitted_curves",
)
DEFAULT_CUTS = ("FR", "IP", "OOP")
LOW_Q_WINDOW = "lamellar_low_q"
TARGET_TEMPERATURES_C = (40.0, 80.0, 150.0)
DETECTION_DELTA_BIC_MIN = 6.0
DETECTION_AREA_SNR_MIN = 3.0
CCL_MIN_POINTS_ACROSS_FWHM = 5.0
CCL_MAX_ABS_STANDARDIZED_RESIDUAL = 5.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine revised nine-frame GIWAXS fits, apply reporting gates, "
            "and make a supervisor-ready QC package."
        )
    )
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "drop40"
            / "analysis_9frame_fast"
        ),
        help="Folder containing final and candidate run directories.",
    )
    parser.add_argument(
        "--input",
        action="append",
        default=[],
        metavar="CUT=PATH",
        help=(
            "Final result directory for a cut; repeat for FR, IP and OOP. "
            "Defaults to CUT_revised under --analysis-root."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results"
            / "drop40"
            / "analysis_9frame_qc_report"
        ),
        help="Destination for combined tables, report, and figures.",
    )
    return parser.parse_args()


def parse_final_inputs(
    specifications: list[str], analysis_root: Path
) -> dict[str, Path]:
    if not specifications:
        return {
            cut: analysis_root / f"{cut}_revised" for cut in DEFAULT_CUTS
        }

    inputs: dict[str, Path] = {}
    for specification in specifications:
        if "=" not in specification:
            raise ValueError(
                f"Invalid --input {specification!r}; expected CUT=PATH."
            )
        cut, path_text = specification.split("=", 1)
        cut = cut.strip().upper()
        if not cut:
            raise ValueError(f"Invalid --input {specification!r}: empty cut.")
        inputs[cut] = Path(path_text).expanduser()
    return inputs


def missing_tables(directory: Path) -> list[str]:
    return [
        f"{name}.csv"
        for name in TABLE_NAMES
        if not (directory / f"{name}.csv").is_file()
    ]


def read_run(directory: Path) -> dict[str, pd.DataFrame]:
    missing = missing_tables(directory)
    if missing:
        raise FileNotFoundError(
            f"{directory} is incomplete; missing {', '.join(missing)}"
        )
    return {
        name: pd.read_csv(directory / f"{name}.csv")
        for name in TABLE_NAMES
    }


def truthy(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )


def numeric(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def add_provenance(
    table: pd.DataFrame, cut: str, source_directory: Path
) -> pd.DataFrame:
    result = table.copy()
    result.insert(0, "cut", cut)
    result.insert(1, "source_run", source_directory.name)
    result.insert(2, "source_directory", str(source_directory.resolve()))
    return result


def bound_tokens(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {
        token.strip()
        for token in re.split(r"[;,]", str(value))
        if token.strip()
    }


def join_quality_on_peaks(
    parameters: pd.DataFrame, quality: pd.DataFrame
) -> pd.DataFrame:
    result = parameters.copy()
    result["_cut_join"] = result["cut"].astype(str)
    result["_source_join"] = result["source_run"].astype(str)
    result["_frame_join"] = result["frame"].astype(str)
    result["_window_join"] = result["window"].astype(str)

    quality_subset = quality.copy()
    quality_subset["_cut_join"] = quality_subset["cut"].astype(str)
    quality_subset["_source_join"] = quality_subset["source_run"].astype(str)
    quality_subset["_frame_join"] = quality_subset["frame"].astype(str)
    quality_subset["_window_join"] = quality_subset["window"].astype(str)
    columns = [
        "_cut_join",
        "_source_join",
        "_frame_join",
        "_window_join",
        "parameters_at_bounds",
        "max_abs_standardized_residual",
        "optimizer_success",
        "quality_flags",
    ]
    available = [column for column in columns if column in quality_subset]
    quality_subset = quality_subset[available].rename(
        columns={
            "quality_flags": "window_quality_flags",
            "optimizer_success": "window_optimizer_success",
        }
    )
    result = result.merge(
        quality_subset,
        how="left",
        on=[
            "_cut_join",
            "_source_join",
            "_frame_join",
            "_window_join",
        ],
        validate="many_to_one",
    )
    return result.drop(
        columns=[
            "_cut_join",
            "_source_join",
            "_frame_join",
            "_window_join",
        ]
    )


def gate_reasons(row: pd.Series, ccl: bool) -> list[str]:
    reasons: list[str] = []
    delta_bic = pd.to_numeric(row.get("delta_bic_peak"), errors="coerce")
    area_snr = pd.to_numeric(row.get("area_snr"), errors="coerce")
    if not np.isfinite(delta_bic) or delta_bic < DETECTION_DELTA_BIC_MIN:
        reasons.append("delta_BIC_peak<6_or_missing")
    if not np.isfinite(area_snr) or area_snr < DETECTION_AREA_SNR_MIN:
        reasons.append("area_SNR<3_or_missing")
    if bool(row.get("center_at_bound", False)):
        reasons.append("centre_at_bound")
    if bool(row.get("fwhm_at_bound", False)):
        reasons.append("FWHM_at_bound")

    if ccl:
        points = pd.to_numeric(
            row.get("points_across_fwhm"), errors="coerce"
        )
        reportable_width = pd.to_numeric(
            row.get("apparent_fwhm_reported_Ainv"), errors="coerce"
        )
        max_residual = pd.to_numeric(
            row.get("window_max_abs_standardized_residual"),
            errors="coerce",
        )
        width_reportable = bool(row.get("_fwhm_reportable_bool", False))
        if (
            not np.isfinite(points)
            or points < CCL_MIN_POINTS_ACROSS_FWHM
        ):
            reasons.append("fewer_than_5_q_steps_across_FWHM")
        if not width_reportable or not np.isfinite(reportable_width):
            reasons.append("fitted_width_not_reportable")
        if (
            not np.isfinite(max_residual)
            or max_residual
            >= CCL_MAX_ABS_STANDARDIZED_RESIDUAL
        ):
            reasons.append("max_abs_standardized_residual>=5_or_missing")
    return list(dict.fromkeys(reasons))


def apply_reporting_gates(
    parameters: pd.DataFrame, quality: pd.DataFrame
) -> pd.DataFrame:
    result = join_quality_on_peaks(parameters, quality)
    result["parameters_at_bounds"] = result.get(
        "parameters_at_bounds", ""
    ).fillna("")

    bound_sets = result["parameters_at_bounds"].map(bound_tokens)
    peak_names = result["peak"].astype(str)
    result["center_at_bound"] = [
        f"{peak}.center" in bounds
        for peak, bounds in zip(peak_names, bound_sets)
    ]
    result["fwhm_at_bound"] = [
        f"{peak}.fwhm" in bounds
        for peak, bounds in zip(peak_names, bound_sets)
    ]
    if "fwhm_reportable" in result:
        result["_fwhm_reportable_bool"] = truthy(
            result["fwhm_reportable"]
        )
    else:
        result["_fwhm_reportable_bool"] = False
    if "optimizer_success" in result:
        result["_optimizer_success_bool"] = truthy(
            result["optimizer_success"]
        )
    elif "window_optimizer_success" in result:
        result["_optimizer_success_bool"] = truthy(
            result["window_optimizer_success"]
        )
    else:
        result["_optimizer_success_bool"] = False

    if "max_abs_standardized_residual" in result:
        result = result.rename(
            columns={
                "max_abs_standardized_residual": (
                    "window_max_abs_standardized_residual"
                )
            }
        )

    detection_reasons = result.apply(
        lambda row: gate_reasons(row, ccl=False), axis=1
    )
    result["detection_gate_pass"] = detection_reasons.map(len).eq(0)
    result["detection_gate_exclusion_reasons"] = detection_reasons.map(
        lambda reasons: ";".join(reasons)
    )

    ccl_reasons = result.apply(
        lambda row: gate_reasons(row, ccl=True), axis=1
    )
    result["apparent_coherence_candidate_gate_pass"] = (
        ccl_reasons.map(len).eq(0)
    )
    result[
        "apparent_coherence_candidate_exclusion_reasons"
    ] = ccl_reasons.map(lambda reasons: ";".join(reasons))

    fwhm = numeric(result["apparent_fwhm_reported_Ainv"])
    points = numeric(result["points_across_fwhm"])
    result["estimated_q_step_Ainv"] = fwhm / points
    result[
        "uncorrected_reciprocal_space_apparent_coherence_length_A"
    ] = np.where(
        result["apparent_coherence_candidate_gate_pass"] & (fwhm > 0),
        2.0 * np.pi / fwhm,
        np.nan,
    )
    result["instrument_correction_status"] = np.where(
        result["apparent_coherence_candidate_gate_pass"],
        "pending; value is not final CCL",
        "not calculated because reporting gate failed",
    )
    return result.drop(
        columns=[
            "_fwhm_reportable_bool",
            "_optimizer_success_bool",
        ],
        errors="ignore",
    )


def selected_rows(model_selection: pd.DataFrame) -> pd.DataFrame:
    if "selected" not in model_selection:
        return model_selection.iloc[0:0].copy()
    return model_selection.loc[truthy(model_selection["selected"])].copy()


def infer_cut(directory: Path) -> str | None:
    match = re.match(r"^(FR|IP|OOP)(?:_|$)", directory.name.upper())
    return match.group(1) if match else None


def candidate_record(
    cut: str,
    directory: Path,
    is_final: bool,
) -> dict[str, object] | None:
    try:
        run = read_run(directory)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None
    model = run["model_selection"]
    selected = selected_rows(model)
    selected = selected.loc[selected["window"].astype(str).eq(LOW_Q_WINDOW)]
    if selected.empty:
        return None
    selected_row = selected.sort_values("bic").iloc[0]

    peaks = run["peak_parameters"]
    peaks = peaks.loc[peaks["window"].astype(str).eq(LOW_Q_WINDOW)]
    curves = run["fitted_curves"]
    curves = curves.loc[curves["window"].astype(str).eq(LOW_Q_WINDOW)]
    q = numeric(curves["q_Ainv"]) if "q_Ainv" in curves else pd.Series()

    return {
        "cut": cut,
        "run": directory.name,
        "source_directory": str(directory.resolve()),
        "is_final_run": is_final,
        "low_q_window": LOW_Q_WINDOW,
        "n_low_q_components": int(peaks["peak"].nunique()),
        "n_frames": int(peaks["frame"].nunique()),
        "selected_profile_family": selected_row.get("profile", ""),
        "eta_shared": selected_row.get("eta_shared", np.nan),
        "selected_model_bic": selected_row.get("bic", np.nan),
        "n_points": selected_row.get("n_points", np.nan),
        "n_parameters": selected_row.get("n_parameters", np.nan),
        "q_min_Ainv": q.min() if not q.empty else np.nan,
        "q_max_Ainv": q.max() if not q.empty else np.nan,
    }


def discover_low_q_candidates(
    analysis_root: Path, final_inputs: dict[str, Path]
) -> pd.DataFrame:
    directories: dict[Path, tuple[str, bool]] = {}
    final_resolved = {
        path.resolve(): cut for cut, path in final_inputs.items()
    }
    for cut, path in final_inputs.items():
        directories[path.resolve()] = (cut, True)
    if analysis_root.is_dir():
        for directory in analysis_root.iterdir():
            if not directory.is_dir():
                continue
            # Background/window sensitivity runs are manual robustness checks,
            # not component-count candidates. They must not silently win an
            # automated BIC comparison against the settled background model.
            if "sensitivity" in directory.name.lower():
                continue
            cut = infer_cut(directory)
            if cut is None:
                continue
            resolved = directory.resolve()
            directories.setdefault(
                resolved, (cut, resolved in final_resolved)
            )

    records = []
    for directory, (cut, is_final) in sorted(
        directories.items(), key=lambda item: (item[1][0], item[0].name)
    ):
        record = candidate_record(cut, directory, is_final)
        if record is not None:
            records.append(record)
    candidates = pd.DataFrame.from_records(records)
    if candidates.empty:
        return candidates

    for column in (
        "selected_model_bic",
        "n_points",
        "q_min_Ainv",
        "q_max_Ainv",
    ):
        candidates[column] = numeric(candidates[column])
    candidates["comparison_signature"] = candidates.apply(
        lambda row: (
            f"n={row['n_points']:.0f};"
            f"q={row['q_min_Ainv']:.6f}-{row['q_max_Ainv']:.6f}"
        ),
        axis=1,
    )

    # A revised run may be a byte-for-byte rerun of the winning candidate.
    # Show that model once, preferring the directory explicitly supplied as
    # the final result.
    candidates["_is_final_sort"] = truthy(
        candidates["is_final_run"]
    ).astype(int)
    candidates["_model_fingerprint"] = candidates.apply(
        lambda row: (
            f"{row['cut']}|{int(row['n_low_q_components'])}|"
            f"{row['selected_profile_family']}|"
            f"{row['selected_model_bic']:.6f}|"
            f"{row['n_points']:.3f}|{row['n_parameters']:.3f}|"
            f"{row['q_min_Ainv']:.6f}|{row['q_max_Ainv']:.6f}"
        ),
        axis=1,
    )
    candidates = (
        candidates.sort_values("_is_final_sort", ascending=False)
        .drop_duplicates("_model_fingerprint", keep="first")
        .drop(columns=["_is_final_sort", "_model_fingerprint"])
        .reset_index(drop=True)
    )
    candidates["comparable_to_final"] = False
    candidates["delta_bic_component_model"] = np.nan
    candidates["component_count_comparison_status"] = (
        "not_compared_no_saved_alternative"
    )
    for cut, cut_rows in candidates.groupby("cut"):
        final_rows = cut_rows.loc[truthy(cut_rows["is_final_run"])]
        if final_rows.empty:
            signature = cut_rows["comparison_signature"].mode().iloc[0]
        else:
            signature = final_rows.iloc[0]["comparison_signature"]
        comparable_index = cut_rows.index[
            cut_rows["comparison_signature"].eq(signature)
        ]
        candidates.loc[comparable_index, "comparable_to_final"] = True
        if len(comparable_index) >= 2:
            best_bic = candidates.loc[
                comparable_index, "selected_model_bic"
            ].min()
            candidates.loc[
                comparable_index, "delta_bic_component_model"
            ] = (
                candidates.loc[comparable_index, "selected_model_bic"]
                - best_bic
            )
            candidates.loc[
                comparable_index, "component_count_comparison_status"
            ] = "compared_by_BIC"
    return candidates.sort_values(
        ["cut", "comparable_to_final", "selected_model_bic"],
        ascending=[True, False, True],
    )


def final_profile_comparison(
    model_selection: pd.DataFrame,
) -> pd.DataFrame:
    low_q = model_selection.loc[
        model_selection["window"].astype(str).eq(LOW_Q_WINDOW)
    ].copy()
    if low_q.empty:
        return low_q
    low_q["bic"] = numeric(low_q["bic"])
    low_q["delta_bic_profile_family"] = (
        low_q["bic"] - low_q.groupby("cut")["bic"].transform("min")
    )
    low_q["best_to_runner_up_delta_bic"] = np.nan
    low_q["profile_family_decision"] = "not_assessable"
    for cut, indices in low_q.groupby("cut").groups.items():
        deltas = (
            low_q.loc[indices, "delta_bic_profile_family"]
            .dropna()
            .sort_values()
            .to_numpy()
        )
        if len(deltas) < 2:
            continue
        margin = float(deltas[1])
        low_q.loc[indices, "best_to_runner_up_delta_bic"] = margin
        low_q.loc[indices, "profile_family_decision"] = (
            "decisive_delta_BIC>=6"
            if margin >= 6.0
            else "not_decisive_delta_BIC<6"
        )
    return low_q


def save_figure(
    figure: plt.Figure, figure_directory: Path, stem: str
) -> None:
    figure.savefig(
        figure_directory / f"{stem}.png",
        dpi=220,
        bbox_inches="tight",
        facecolor="white",
    )
    figure.savefig(
        figure_directory / f"{stem}.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)


def plot_low_q_decisions(
    candidates: pd.DataFrame,
    profile_comparison: pd.DataFrame,
    figure_directory: Path,
) -> bool:
    cuts = sorted(
        set(candidates.get("cut", pd.Series(dtype=str)).dropna())
        | set(profile_comparison.get("cut", pd.Series(dtype=str)).dropna())
    )
    if not cuts:
        return False
    figure, axes = plt.subplots(
        len(cuts),
        2,
        figsize=(12.5, max(3.8, 3.8 * len(cuts))),
        squeeze=False,
        constrained_layout=True,
    )
    family_labels = {
        "gaussian": "Gaussian",
        "lorentzian": "Lorentzian",
        "pseudo_voigt": "pseudo-Voigt",
    }
    for row_index, cut in enumerate(cuts):
        candidate_axis, profile_axis = axes[row_index]
        cut_candidates = candidates.loc[
            candidates["cut"].eq(cut)
            & truthy(candidates["comparable_to_final"])
        ].copy()
        if cut_candidates.empty:
            candidate_axis.text(
                0.5,
                0.5,
                "No comparable component-count runs",
                ha="center",
                va="center",
                transform=candidate_axis.transAxes,
            )
            candidate_axis.set_axis_off()
        elif len(cut_candidates) == 1:
            only = cut_candidates.iloc[0]
            component_word = (
                "component"
                if int(only["n_low_q_components"]) == 1
                else "components"
            )
            candidate_axis.text(
                0.5,
                0.57,
                (
                    f"Final: {int(only['n_low_q_components'])} "
                    f"{component_word}\n"
                    f"BIC = {float(only['selected_model_bic']):.1f}"
                ),
                ha="center",
                va="center",
                fontsize=14,
                fontweight="bold",
                transform=candidate_axis.transAxes,
            )
            candidate_axis.text(
                0.5,
                0.30,
                (
                    "No saved alternative component-count run\n"
                    "in this report package"
                ),
                ha="center",
                va="center",
                fontsize=10,
                color="#555555",
                transform=candidate_axis.transAxes,
            )
            candidate_axis.set_title(
                f"{cut}: low-q component-count record"
            )
            candidate_axis.set_axis_off()
        else:
            cut_candidates = cut_candidates.sort_values(
                "n_low_q_components"
            )
            labels = [
                (
                    f"{int(row.n_low_q_components)} comp."
                    + ("\nfinal" if bool(row.is_final_run) else "")
                )
                for row in cut_candidates.itertuples()
            ]
            values = numeric(
                cut_candidates["delta_bic_component_model"]
            ).to_numpy()
            colors = [
                "#176B87" if bool(value) else "#A9B8C6"
                for value in truthy(cut_candidates["is_final_run"])
            ]
            bars = candidate_axis.bar(
                np.arange(len(values)), values, color=colors
            )
            for bar, value, bic in zip(
                bars,
                values,
                numeric(cut_candidates["selected_model_bic"]),
            ):
                candidate_axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    f"Δ={value:.1f}\nBIC={bic:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
            candidate_axis.set_xticks(np.arange(len(labels)), labels)
            candidate_axis.tick_params(axis="x", labelsize=8)
            candidate_axis.axhline(
                6,
                color="#B3413E",
                linestyle="--",
                linewidth=1,
                label="ΔBIC = 6",
            )
            candidate_axis.set_ylabel("ΔBIC (component model)")
            candidate_axis.set_title(
                f"{cut}: low-q component-count decision"
            )
            candidate_axis.legend(frameon=False, fontsize=8)
            candidate_axis.spines[["top", "right"]].set_visible(False)

        profiles = profile_comparison.loc[
            profile_comparison["cut"].eq(cut)
        ].copy()
        if profiles.empty:
            profile_axis.text(
                0.5,
                0.5,
                "No final profile-family table",
                ha="center",
                va="center",
                transform=profile_axis.transAxes,
            )
            profile_axis.set_axis_off()
        else:
            profiles = profiles.sort_values(
                "delta_bic_profile_family"
            )
            labels = [
                family_labels.get(str(profile), str(profile))
                for profile in profiles["profile"]
            ]
            values = numeric(
                profiles["delta_bic_profile_family"]
            ).to_numpy()
            selected = (
                truthy(profiles["selected"])
                if "selected" in profiles
                else pd.Series(False, index=profiles.index)
            )
            colors = [
                "#D97841" if bool(value) else "#D5D9DD"
                for value in selected
            ]
            bars = profile_axis.bar(
                np.arange(len(values)), values, color=colors
            )
            for bar, value in zip(bars, values):
                profile_axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    f"{value:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
            profile_axis.set_xticks(np.arange(len(labels)), labels)
            profile_axis.tick_params(axis="x", labelsize=8)
            profile_axis.axhline(
                6,
                color="#B3413E",
                linestyle="--",
                linewidth=1,
            )
            profile_axis.set_ylabel("ΔBIC (profile family)")
            profile_axis.set_title(
                f"{cut}: one family selected across all 9 frames"
            )
            profile_axis.spines[["top", "right"]].set_visible(False)

    figure.suptitle(
        "Low-q model decisions (ΔBIC = 0 is preferred)",
        fontsize=14,
        fontweight="bold",
    )
    save_figure(figure, figure_directory, "low_q_model_decision")
    return True


def component_columns(curves: pd.DataFrame) -> list[str]:
    columns = []
    for column in curves.columns:
        if not column.startswith("component_"):
            continue
        values = numeric(curves[column])
        if values.notna().any() and np.nansum(np.abs(values)) > 0:
            columns.append(column)
    return columns


def closest_temperature(
    available: Iterable[float], target: float
) -> float | None:
    values = np.asarray(list(available), dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None
    return float(values[np.argmin(np.abs(values - target))])


def plot_low_q_overlays(
    curves: pd.DataFrame, figure_directory: Path
) -> bool:
    low_q = curves.loc[curves["window"].astype(str).eq(LOW_Q_WINDOW)].copy()
    if low_q.empty:
        return False
    cuts = [cut for cut in DEFAULT_CUTS if cut in set(low_q["cut"])]
    cuts += sorted(set(low_q["cut"]) - set(cuts))
    if not cuts:
        return False

    figure = plt.figure(
        figsize=(15.2, 3.6 * len(cuts)),
        constrained_layout=True,
    )
    grid = figure.add_gridspec(
        nrows=2 * len(cuts),
        ncols=len(TARGET_TEMPERATURES_C),
        height_ratios=[
            value
            for _ in cuts
            for value in (3.0, 1.0)
        ],
    )
    fit_colours = {
        "FR": "#176B87",
        "IP": "#7A5195",
        "OOP": "#D97841",
    }
    component_palette = [
        "#2A9D8F",
        "#E9C46A",
        "#E76F51",
        "#6A4C93",
        "#8AB17D",
    ]

    for cut_index, cut in enumerate(cuts):
        cut_data = low_q.loc[low_q["cut"].eq(cut)].copy()
        available_temperatures = numeric(cut_data["temperature_C"]).unique()
        for column_index, target in enumerate(TARGET_TEMPERATURES_C):
            main_axis = figure.add_subplot(
                grid[2 * cut_index, column_index]
            )
            residual_axis = figure.add_subplot(
                grid[2 * cut_index + 1, column_index],
                sharex=main_axis,
            )
            actual = closest_temperature(available_temperatures, target)
            if actual is None:
                main_axis.text(
                    0.5,
                    0.5,
                    f"No {target:g} °C frame",
                    ha="center",
                    va="center",
                    transform=main_axis.transAxes,
                )
                main_axis.set_axis_off()
                residual_axis.set_axis_off()
                continue
            panel = cut_data.loc[
                np.isclose(
                    numeric(cut_data["temperature_C"]),
                    actual,
                    atol=1e-6,
                )
            ].sort_values("q_Ainv")
            if panel.empty:
                main_axis.set_axis_off()
                residual_axis.set_axis_off()
                continue

            q = numeric(panel["q_Ainv"]).to_numpy()
            observed = numeric(panel["observed_intensity"]).to_numpy()
            total = numeric(panel["total_fit_intensity"]).to_numpy()
            background = numeric(
                panel["background_intensity"]
            ).to_numpy()
            main_axis.plot(
                q,
                observed,
                "o",
                markersize=2.8,
                color="#252525",
                alpha=0.72,
                label="Raw data",
                zorder=3,
            )
            main_axis.plot(
                q,
                total,
                color=fit_colours.get(cut, "#176B87"),
                linewidth=2.2,
                label="Total fit",
                zorder=4,
            )
            main_axis.plot(
                q,
                background,
                color="#777777",
                linestyle="--",
                linewidth=1.3,
                label="Background",
                zorder=2,
            )
            for component_index, component in enumerate(
                component_columns(panel)
            ):
                values = numeric(panel[component]).to_numpy()
                label = component.removeprefix("component_")
                main_axis.plot(
                    q,
                    background + values,
                    linewidth=1.1,
                    color=component_palette[
                        component_index % len(component_palette)
                    ],
                    alpha=0.9,
                    label=label,
                    zorder=1,
                )

            frame = panel["frame"].astype(str).iloc[0]
            main_axis.set_title(
                f"{cut} — {target:g} °C nominal (scan {frame})",
                fontsize=10.5,
            )
            if column_index == 0:
                main_axis.set_ylabel("Normalised intensity (a.u.)")
                residual_axis.set_ylabel("Std.\nresidual")
            main_axis.tick_params(axis="x", labelbottom=False)
            main_axis.spines[["top", "right"]].set_visible(False)

            standardized = numeric(
                panel["standardized_residual"]
            ).to_numpy()
            residual_axis.axhspan(
                -3, 3, color="#B8D8D8", alpha=0.25, zorder=0
            )
            residual_axis.axhline(0, color="#555555", linewidth=0.8)
            residual_axis.axhline(
                5, color="#B3413E", linestyle=":", linewidth=0.8
            )
            residual_axis.axhline(
                -5, color="#B3413E", linestyle=":", linewidth=0.8
            )
            residual_axis.plot(
                q,
                standardized,
                "o-",
                color="#555555",
                markersize=2.2,
                linewidth=0.8,
            )
            finite = standardized[np.isfinite(standardized)]
            extent = max(
                5.5,
                min(12.0, np.nanmax(np.abs(finite)) * 1.15)
                if finite.size
                else 5.5,
            )
            residual_axis.set_ylim(-extent, extent)
            residual_axis.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
            residual_axis.spines[["top", "right"]].set_visible(False)
            if column_index == 0:
                main_axis.legend(
                    frameon=False,
                    fontsize=7.2,
                    ncol=2,
                    loc="best",
                )

    figure.suptitle(
        "Low-q fit diagnostics: raw data, total model, components and residuals",
        fontsize=14,
        fontweight="bold",
    )
    save_figure(figure, figure_directory, "low_q_overlays_40_80_150")
    return True


def build_gate_summary(parameters: pd.DataFrame) -> pd.DataFrame:
    grouped = parameters.groupby(
        ["cut", "window", "window_title"], dropna=False
    )
    summary = grouped.agg(
        fitted_peak_rows=("peak", "size"),
        detection_gate_passes=("detection_gate_pass", "sum"),
        apparent_coherence_candidates=(
            "apparent_coherence_candidate_gate_pass",
            "sum",
        ),
        unique_peaks=("peak", "nunique"),
        frames=("frame", "nunique"),
    ).reset_index()
    summary["detection_gate_pass_fraction"] = (
        summary["detection_gate_passes"]
        / summary["fitted_peak_rows"]
    )
    summary["apparent_coherence_candidate_fraction"] = (
        summary["apparent_coherence_candidates"]
        / summary["fitted_peak_rows"]
    )
    return summary


def plot_gate_summary(
    summary: pd.DataFrame, figure_directory: Path
) -> bool:
    if summary.empty:
        return False
    labels = (
        summary["cut"].astype(str)
        + " · "
        + summary["window"].astype(str)
    )
    values = summary[
        [
            "detection_gate_pass_fraction",
            "apparent_coherence_candidate_fraction",
        ]
    ].to_numpy(dtype=float)
    figure_height = max(4.0, 0.34 * len(summary) + 1.8)
    figure, axis = plt.subplots(
        figsize=(9.0, figure_height), constrained_layout=True
    )
    image = axis.imshow(
        values,
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap="Blues",
    )
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            axis.text(
                column,
                row,
                (
                    f"{values[row, column] * 100:.0f}%\n"
                    f"({int(summary.iloc[row][['detection_gate_passes', 'apparent_coherence_candidates'][column]])}"
                    f"/{int(summary.iloc[row]['fitted_peak_rows'])})"
                ),
                ha="center",
                va="center",
                color="white" if values[row, column] > 0.55 else "#202020",
                fontsize=8,
            )
    axis.set_xticks(
        [0, 1],
        [
            "Detection claim gate",
            "Apparent-coherence candidate gate",
        ],
    )
    axis.set_yticks(np.arange(len(labels)), labels)
    axis.set_title(
        "Automated reporting-gate pass rates\n"
        "(manual background/window review still required)",
        fontweight="bold",
    )
    colourbar = figure.colorbar(image, ax=axis, shrink=0.7)
    colourbar.set_label("Fraction passing")
    save_figure(figure, figure_directory, "reporting_gate_summary")
    return True


def markdown_table(
    table: pd.DataFrame, columns: list[str], max_rows: int | None = None
) -> str:
    available = [column for column in columns if column in table]
    if not available or table.empty:
        return "_None._"
    subset = table[available]
    if max_rows is not None:
        subset = subset.head(max_rows)

    def display(value: object) -> str:
        if pd.isna(value):
            return "—"
        if isinstance(value, (float, np.floating)):
            if abs(value) >= 1000:
                return f"{value:.1f}"
            return f"{value:.4g}"
        if isinstance(value, (bool, np.bool_)):
            return "yes" if value else "no"
        return str(value).replace("|", r"\|").replace("\n", " ")

    header = "| " + " | ".join(available) + " |"
    separator = "| " + " | ".join("---" for _ in available) + " |"
    rows = [
        "| "
        + " | ".join(display(value) for value in row)
        + " |"
        for row in subset.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def candidate_peak_summary(
    candidates: pd.DataFrame,
) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    grouped = candidates.groupby(
        ["cut", "window", "peak", "peak_label"], dropna=False
    )
    summary = grouped.agg(
        passing_frames=("frame", "nunique"),
        nominal_temperature_min_C=("temperature_C", "min"),
        nominal_temperature_max_C=("temperature_C", "max"),
        q0_min_Ainv=("q0_fit_Ainv", "min"),
        q0_max_Ainv=("q0_fit_Ainv", "max"),
        apparent_fwhm_min_Ainv=("apparent_fwhm_reported_Ainv", "min"),
        apparent_fwhm_max_Ainv=("apparent_fwhm_reported_Ainv", "max"),
        uncorrected_apparent_length_min_A=(
            "uncorrected_reciprocal_space_apparent_coherence_length_A",
            "min",
        ),
        uncorrected_apparent_length_max_A=(
            "uncorrected_reciprocal_space_apparent_coherence_length_A",
            "max",
        ),
    ).reset_index()
    return summary.sort_values(["cut", "window", "peak"])


def build_report(
    final_inputs: dict[str, Path],
    model_selection: pd.DataFrame,
    quality: pd.DataFrame,
    parameters: pd.DataFrame,
    candidates: pd.DataFrame,
    profile_comparison: pd.DataFrame,
    gate_summary: pd.DataFrame,
    figure_status: dict[str, bool],
) -> str:
    passing = parameters.loc[
        truthy(parameters["apparent_coherence_candidate_gate_pass"])
    ].copy()
    candidate_summary = candidate_peak_summary(passing)

    quality_copy = quality.copy()
    quality_copy["_optimizer"] = (
        truthy(quality_copy["optimizer_success"])
        if "optimizer_success" in quality_copy
        else False
    )
    max_residual = numeric(
        quality_copy.get(
            "max_abs_standardized_residual",
            pd.Series(np.nan, index=quality_copy.index),
        )
    )
    bound_windows = int(
        quality_copy.get(
            "parameters_at_bounds",
            pd.Series("", index=quality_copy.index),
        )
        .fillna("")
        .astype(str)
        .str.strip()
        .ne("")
        .sum()
    )
    large_residual_windows = int(
        (max_residual >= CCL_MAX_ABS_STANDARDIZED_RESIDUAL).sum()
    )
    failed_optimizers = int((~quality_copy["_optimizer"]).sum())

    bootstrap_requested = numeric(
        parameters.get(
            "bootstrap_requested",
            pd.Series(0, index=parameters.index),
        )
    )
    maximum_bootstrap = (
        int(bootstrap_requested.max())
        if bootstrap_requested.notna().any()
        else 0
    )
    if maximum_bootstrap == 0:
        bootstrap_note = (
            "These are fast, no-bootstrap screening results. Run the settled "
            "models with 100 moving-block bootstrap repeats before using "
            "uncertainty intervals in the dissertation."
        )
    else:
        bootstrap_note = (
            f"The source tables request up to {maximum_bootstrap} bootstrap "
            "repeats per peak/frame. Check `bootstrap_successes` before "
            "using each interval."
        )

    profile_selected = selected_rows(profile_comparison)
    profile_columns = [
        "cut",
        "profile",
        "eta_shared",
        "bic",
        "best_to_runner_up_delta_bic",
        "profile_family_decision",
    ]
    component_columns_report = [
        "cut",
        "run",
        "is_final_run",
        "n_low_q_components",
        "selected_profile_family",
        "selected_model_bic",
        "delta_bic_component_model",
        "comparable_to_final",
        "component_count_comparison_status",
    ]
    gate_columns = [
        "cut",
        "window",
        "fitted_peak_rows",
        "detection_gate_passes",
        "apparent_coherence_candidates",
    ]
    apparent_columns = [
        "cut",
        "window",
        "peak",
        "passing_frames",
        "nominal_temperature_min_C",
        "nominal_temperature_max_C",
        "apparent_fwhm_min_Ainv",
        "apparent_fwhm_max_Ainv",
        "uncorrected_apparent_length_min_A",
        "uncorrected_apparent_length_max_A",
    ]

    input_lines = "\n".join(
        f"- **{cut}:** `{path.resolve()}`"
        for cut, path in final_inputs.items()
    )
    figure_lines = []
    if figure_status.get("low_q_model_decision"):
        figure_lines.append(
            "- [Low-q model decision](figures/low_q_model_decision.png)"
        )
    if figure_status.get("low_q_overlays"):
        figure_lines.append(
            "- [40/80/150 °C low-q overlays]"
            "(figures/low_q_overlays_40_80_150.png)"
        )
    if figure_status.get("reporting_gate_summary"):
        figure_lines.append(
            "- [Reporting-gate summary]"
            "(figures/reporting_gate_summary.png)"
        )
    figures_markdown = "\n".join(figure_lines) or "_No figures generated._"

    lines = [
        "# Nine-frame GIWAXS peak-fit QC report",
        "",
        (
            "**Temperatures in this report are nominal temperatures.** "
            "The series is 40, 50, 55, 60, 65, 70, 75, 80 and 150 °C."
        ),
        "",
        (
            "**Scientific boundary:** no value in this report is a final "
            "crystalline coherence length (CCL). Where the automated gates "
            "pass, `2π/FWHM` is shown only as an **uncorrected "
            "reciprocal-space apparent coherence length**. Instrumental "
            "broadening correction, a defensible peak assignment, and a "
            "manual background/window sensitivity check are still required."
        ),
        "",
        "## Final source runs",
        "",
        input_lines,
        "",
        "## Supervisor-ready figures",
        "",
        figures_markdown,
        "",
        "## Low-q model decision",
        "",
        (
            "The component-count comparison uses only candidate runs with "
            "the same number of fitted points and the same low-q range as "
            "the final run. If only one component model is listed, the "
            "package records the final choice but does not claim to have "
            "retested its component count. The profile-family comparison "
            "fits Gaussian, "
            "Lorentzian and pseudo-Voigt families across all nine frames; "
            "it is not a majority vote over frames. Lower BIC is preferred, "
            "and ΔBIC ≥ 6 is treated as substantial evidence. A selected "
            "family marked `not_decisive_delta_BIC<6` should be retained for "
            "consistent tracking, but interpreted as practically tied with "
            "the runner-up rather than uniquely established."
        ),
        "",
        "### Comparable component models",
        "",
        markdown_table(candidates, component_columns_report),
        "",
        "### Profile family selected for each final cut",
        "",
        markdown_table(profile_selected, profile_columns),
        "",
        "## Automated reporting gates",
        "",
        (
            "A per-frame peak may be claimed as detected only when "
            "`delta_BIC_peak ≥ 6`, `area_SNR ≥ 3`, and neither that peak's "
            "centre nor its FWHM is at a bound."
        ),
        "",
        (
            "An apparent-coherence candidate must also have at least five "
            "q steps across its FWHM, a reportable fitted width, no FWHM "
            "bound, and a window maximum absolute standardised residual "
            "below 5. Passing this automated gate does **not** test "
            "background/window sensitivity."
        ),
        "",
        markdown_table(gate_summary, gate_columns),
        "",
        "## Apparent-coherence candidates",
        "",
        (
            "The following summary includes only frames that pass all "
            "automated gates. The underlying per-frame rows are in "
            "`ccl_candidates.csv`."
        ),
        "",
        markdown_table(candidate_summary, apparent_columns),
        "",
        "## Fit-level warnings",
        "",
        f"- Optimiser failures: **{failed_optimizers}** fit windows.",
        f"- Fit windows with at least one parameter at a bound: **{bound_windows}**.",
        (
            "- Fit windows with maximum absolute standardised residual "
            f"≥ 5: **{large_residual_windows}**."
        ),
        f"- {bootstrap_note}",
        "",
        "## Required manual review before dissertation use",
        "",
        (
            "1. Inspect every retained overlay and residual plot, especially "
            "40, 80 and 150 °C."
        ),
        (
            "2. Repeat each proposed width with sensible alternative "
            "background orders and slightly shifted fit windows. Reject a "
            "width if the result is not stable."
        ),
        (
            "3. Confirm the physical peak assignment and that the component "
            "is not merely a mathematical shoulder."
        ),
        (
            "4. Correct acceptable sample widths for instrumental broadening "
            "using a resolvable LaB6 line measured and integrated in matching "
            "geometry. Do not subtract FWHM values directly."
        ),
        (
            "5. Only after the model is frozen, run 100 moving-block bootstrap "
            "repeats per frame and report uncertainty intervals."
        ),
        "",
        "## Output tables",
        "",
        "- `model_selection_combined.csv`",
        "- `fit_quality_combined.csv`",
        "- `peak_parameters_combined_with_gates.csv`",
        "- `fitted_curves_combined.csv`",
        "- `low_q_model_comparison.csv`",
        "- `low_q_profile_family_comparison.csv`",
        "- `reporting_gate_summary.csv`",
        "- `detection_claim_candidates.csv`",
        "- `ccl_candidates.csv` (candidate name retained for workflow "
        "compatibility; values inside are explicitly non-final)",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    analysis_root = args.analysis_root.expanduser()
    final_inputs = parse_final_inputs(args.input, analysis_root)
    output = args.output.expanduser()
    figure_directory = output / "figures"
    output.mkdir(parents=True, exist_ok=True)
    figure_directory.mkdir(parents=True, exist_ok=True)

    loaded: dict[str, dict[str, pd.DataFrame]] = {}
    incomplete: list[str] = []
    for cut, directory in final_inputs.items():
        try:
            loaded[cut] = read_run(directory)
        except FileNotFoundError as exc:
            incomplete.append(str(exc))
    if incomplete:
        raise SystemExit(
            "Final result inputs are incomplete:\n- "
            + "\n- ".join(incomplete)
        )

    combined: dict[str, pd.DataFrame] = {}
    for table_name in TABLE_NAMES:
        combined[table_name] = pd.concat(
            [
                add_provenance(
                    loaded[cut][table_name], cut, final_inputs[cut]
                )
                for cut in final_inputs
            ],
            ignore_index=True,
            sort=False,
        )

    gated_parameters = apply_reporting_gates(
        combined["peak_parameters"], combined["fit_quality"]
    )
    candidates = discover_low_q_candidates(analysis_root, final_inputs)
    profile_comparison = final_profile_comparison(
        combined["model_selection"]
    )
    gate_summary = build_gate_summary(gated_parameters)

    combined["model_selection"].to_csv(
        output / "model_selection_combined.csv", index=False
    )
    combined["fit_quality"].to_csv(
        output / "fit_quality_combined.csv", index=False
    )
    gated_parameters.to_csv(
        output / "peak_parameters_combined_with_gates.csv", index=False
    )
    combined["fitted_curves"].to_csv(
        output / "fitted_curves_combined.csv", index=False
    )
    candidates.to_csv(output / "low_q_model_comparison.csv", index=False)
    profile_comparison.to_csv(
        output / "low_q_profile_family_comparison.csv", index=False
    )
    gate_summary.to_csv(
        output / "reporting_gate_summary.csv", index=False
    )
    gated_parameters.loc[
        truthy(gated_parameters["detection_gate_pass"])
    ].to_csv(output / "detection_claim_candidates.csv", index=False)
    gated_parameters.loc[
        truthy(
            gated_parameters[
                "apparent_coherence_candidate_gate_pass"
            ]
        )
    ].to_csv(output / "ccl_candidates.csv", index=False)

    figure_status = {
        "low_q_model_decision": plot_low_q_decisions(
            candidates, profile_comparison, figure_directory
        ),
        "low_q_overlays": plot_low_q_overlays(
            combined["fitted_curves"], figure_directory
        ),
        "reporting_gate_summary": plot_gate_summary(
            gate_summary, figure_directory
        ),
    }
    report = build_report(
        final_inputs=final_inputs,
        model_selection=combined["model_selection"],
        quality=combined["fit_quality"],
        parameters=gated_parameters,
        candidates=candidates,
        profile_comparison=profile_comparison,
        gate_summary=gate_summary,
        figure_status=figure_status,
    )
    (output / "QC_REPORT.md").write_text(report, encoding="utf-8")

    print(f"QC package written to: {output.resolve()}")
    print(
        "No final CCL values were produced; any 2π/FWHM values are "
        "explicitly labelled uncorrected apparent coherence lengths."
    )


if __name__ == "__main__":
    main()
