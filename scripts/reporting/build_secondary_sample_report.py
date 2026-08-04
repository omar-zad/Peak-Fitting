#!/usr/bin/env python3
"""Build conservative, dissertation-facing tables for a secondary scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


CUTS = ("FR", "IP", "OOP")
ACCEPTED_FREE_STATUSES = {"strong", "detected"}
ACCEPTED_AREA_STATUSES = ACCEPTED_FREE_STATUSES | {"detected_fixed_shape"}


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def bound_contains(value: object, parameter: str) -> bool:
    if pd.isna(value):
        return False
    names = {item.strip() for item in str(value).split(";") if item.strip()}
    return parameter in names


def exclusion_reasons(row: pd.Series, *, for_position: bool) -> str:
    reasons: list[str] = []
    accepted_statuses = ACCEPTED_FREE_STATUSES if for_position else ACCEPTED_AREA_STATUSES
    if row["detection_status"] not in accepted_statuses:
        reasons.append(f"status={row['detection_status']}")
    if not np.isfinite(row["delta_bic_peak"]) or row["delta_bic_peak"] < 6:
        reasons.append("delta_BIC<6_or_missing")
    if not np.isfinite(row["area_snr"]) or row["area_snr"] < 3:
        reasons.append("area_SNR<3_or_missing")
    if not bool(row["optimizer_success"]):
        reasons.append("optimizer_failed")
    if bool(row["center_at_bound"]):
        reasons.append("centre_at_bound")
    if bool(row.get("reference_like_component", False)):
        reasons.append("reference_like_component")
    if bool(row.get("position_only", False)) and not for_position:
        reasons.append("position_only_no_area_claim")
    if bool(row["fwhm_at_bound"]) and not (
        for_position and bool(row.get("position_only", False))
    ):
        reasons.append("FWHM_at_bound")
    if for_position:
        if not bool(row["position_reportable"]):
            reasons.append("position_not_reportable")
        if not np.isfinite(row["q0_reported_Ainv"]):
            reasons.append("q0_reported_missing")
    return ";".join(dict.fromkeys(reasons))


def load_cut(results_root: Path, cut: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    cut_dir = results_root / cut
    parameter_path = cut_dir / "peak_parameters.csv"
    quality_path = cut_dir / "fit_quality.csv"
    selection_path = cut_dir / "model_selection.csv"
    for path in (parameter_path, quality_path, selection_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    parameters = pd.read_csv(parameter_path)
    quality = pd.read_csv(quality_path)
    selection = pd.read_csv(selection_path)
    parameters["cut"] = cut
    quality_columns = [
        "frame",
        "window",
        "optimizer_success",
        "parameters_at_bounds",
        "quality_flags",
        "max_abs_standardized_residual",
        "reduced_chi2_empirical",
    ]
    missing_quality = sorted(set(quality_columns) - set(quality.columns))
    if missing_quality:
        raise ValueError(
            f"{quality_path} is missing: {', '.join(missing_quality)}"
        )
    quality = quality.loc[:, quality_columns].rename(
        columns={"quality_flags": "window_quality_flags"}
    )
    combined = parameters.merge(
        quality,
        on=["frame", "window"],
        how="left",
        validate="many_to_one",
    )
    combined["optimizer_success"] = truthy(combined["optimizer_success"])
    combined["position_reportable"] = truthy(combined["position_reportable"])
    for optional_flag in ("position_only", "reference_like_component"):
        if optional_flag not in combined.columns:
            combined[optional_flag] = False
        else:
            combined[optional_flag] = truthy(combined[optional_flag])
    for column in (
        "delta_bic_peak",
        "area_snr",
        "q0_reported_Ainv",
        "q0_ci95_low_Ainv",
        "q0_ci95_high_Ainv",
        "d_reported_A",
        "area_total",
        "area_total_ci95_low",
        "area_total_ci95_high",
        "area_in_window",
        "area_in_window_ci95_low",
        "area_in_window_ci95_high",
    ):
        combined[column] = pd.to_numeric(combined[column], errors="coerce")

    combined["center_at_bound"] = combined.apply(
        lambda row: bound_contains(
            row["parameters_at_bounds"], f"{row['peak']}.center"
        ),
        axis=1,
    )
    combined["fwhm_at_bound"] = combined.apply(
        lambda row: bound_contains(
            row["parameters_at_bounds"], f"{row['peak']}.fwhm"
        ),
        axis=1,
    )
    combined["detection_gate_pass"] = (
        combined["detection_status"].isin(ACCEPTED_AREA_STATUSES)
        & (combined["delta_bic_peak"] >= 6)
        & (combined["area_snr"] >= 3)
        & combined["optimizer_success"]
        & ~combined["center_at_bound"]
        & ~combined["fwhm_at_bound"]
        & ~combined["position_only"]
        & ~combined["reference_like_component"]
    )
    combined["position_claim_pass"] = (
        combined["detection_status"].isin(ACCEPTED_FREE_STATUSES)
        & (combined["delta_bic_peak"] >= 6)
        & (combined["area_snr"] >= 3)
        & combined["optimizer_success"]
        & ~combined["center_at_bound"]
        & (~combined["fwhm_at_bound"] | combined["position_only"])
        & ~combined["reference_like_component"]
        & combined["position_reportable"]
        & combined["q0_reported_Ainv"].notna()
    )
    combined["detection_exclusion_reasons"] = combined.apply(
        lambda row: exclusion_reasons(row, for_position=False), axis=1
    )
    combined["position_exclusion_reasons"] = combined.apply(
        lambda row: exclusion_reasons(row, for_position=True), axis=1
    )

    selection["cut"] = cut
    selection["selected"] = truthy(selection["selected"])
    return combined, selection


def write_report(
    path: Path,
    scan: str,
    sample_label: str,
    combined: pd.DataFrame,
    accepted_positions: pd.DataFrame,
    accepted_areas: pd.DataFrame,
    tentative: pd.DataFrame,
) -> None:
    lines = [
        f"# Secondary-sample peak-fit report: {sample_label}",
        "",
        f"Scan: **{scan}**",
        "",
        "## Scientific scope",
        "",
        "This report is for peak presence, fitted q/d position and integrated "
        "area. Apparent FWHM values are retained only as fit diagnostics; no "
        "CCL, crystallite-size or strain result is claimed.",
        "",
        "An accepted position requires a freely resolved `strong` or `detected` "
        "component, ΔBIC ≥ 6, area SNR ≥ 3, optimiser success, no centre "
        "bound, `position_reportable = True`, and a nonblank "
        "`q0_reported_Ainv`. An explicitly configured position-only line may "
        "retain q when only its under-resolved FWHM is at the resolution floor; "
        "its width and area remain excluded.",
        "",
        "A fixed-shape component may support an area/presence claim when its "
        "detection gate passes, but it cannot provide an independently measured "
        "position.",
        "Reference-like shoulder components are excluded from both assignment "
        "and accepted result tables even when they improve the neighbouring "
        "sharp-peak fit.",
        "",
        "## Counts",
        "",
        f"- Fitted component rows: **{len(combined)}**",
        f"- Accepted peak-position rows: **{len(accepted_positions)}**",
        f"- Accepted peak-area/presence rows: **{len(accepted_areas)}**",
        f"- Tentative or excluded rows: **{len(tentative)}**",
        "",
        "## Required visual review",
        "",
        "Before using a row in the dissertation, inspect its fit overlay and "
        "residual panel. Reject mathematical shoulders, background-sensitive "
        "components, peak swapping, and small apparent movements that are not "
        "supported by their confidence intervals and adjacent frames.",
        "",
        "## Output tables",
        "",
        "- `accepted_peak_positions.csv`: values suitable for a conservative q/d table.",
        "- `accepted_peak_areas.csv`: accepted presence/area rows, including clearly labelled fixed-shape detections.",
        "- `tentative_or_excluded_peaks.csv`: rows that must remain tentative or unreported.",
        "- `all_peak_fits_with_gates.csv`: complete audit trail, including numerical fits that were not accepted.",
        "- `selected_profile_families.csv`: profile family selected or locked for each cut/window.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine FR/IP/OOP fit outputs with conservative reporting gates."
    )
    parser.add_argument("--scan", required=True)
    parser.add_argument("--sample-label", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    manifest_path = args.manifest.expanduser().resolve()
    results_root = args.results_root.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(manifest_path)
    manifest_columns = [
        column
        for column in (
            "frame",
            "acquisition_time_s",
            "measured_temperature_C",
            "phase",
            "time_at_high_temperature_s",
            "key_frame_role",
        )
        if column in manifest.columns
    ]
    manifest = manifest.loc[:, manifest_columns].copy()

    all_parameters = []
    all_selection = []
    for cut in CUTS:
        parameters, selection = load_cut(results_root, cut)
        all_parameters.append(parameters)
        all_selection.append(selection)
    combined = pd.concat(all_parameters, ignore_index=True)
    selection = pd.concat(all_selection, ignore_index=True)
    combined = combined.drop(columns=["temperature_C"], errors="ignore").merge(
        manifest,
        on="frame",
        how="left",
        validate="many_to_one",
    )
    combined.insert(0, "scan", int(args.scan))
    combined.insert(1, "sample_label", args.sample_label)
    for optional_column in (
        "acquisition_time_s",
        "measured_temperature_C",
        "phase",
        "time_at_high_temperature_s",
        "key_frame_role",
    ):
        if optional_column not in combined.columns:
            combined[optional_column] = np.nan

    accepted_positions = combined.loc[combined["position_claim_pass"]].copy()
    accepted_areas = combined.loc[combined["detection_gate_pass"]].copy()
    tentative = combined.loc[~combined["position_claim_pass"]].copy()
    sort_columns = ["cut", "frame", "window", "peak"]
    for table in (combined, accepted_positions, accepted_areas, tentative):
        table.sort_values(sort_columns, inplace=True)

    position_columns = [
        "scan",
        "sample_label",
        "cut",
        "frame",
        "acquisition_time_s",
        "measured_temperature_C",
        "phase",
        "key_frame_role",
        "window",
        "peak",
        "peak_label",
        "position_only",
        "reference_like_component",
        "profile",
        "eta_shared",
        "detection_status",
        "delta_bic_peak",
        "area_snr",
        "q0_reported_Ainv",
        "q0_ci95_low_Ainv",
        "q0_ci95_high_Ainv",
        "d_reported_A",
        "area_in_window",
        "area_in_window_ci95_low",
        "area_in_window_ci95_high",
    ]
    area_columns = position_columns[:-7] + [
        "position_reportable",
        "q0_reported_Ainv",
        "area_total",
        "area_total_ci95_low",
        "area_total_ci95_high",
        "area_in_window",
        "area_in_window_ci95_low",
        "area_in_window_ci95_high",
    ]
    accepted_positions.loc[:, position_columns].to_csv(
        outdir / "accepted_peak_positions.csv", index=False
    )
    accepted_areas.loc[:, area_columns].to_csv(
        outdir / "accepted_peak_areas.csv", index=False
    )
    tentative.to_csv(outdir / "tentative_or_excluded_peaks.csv", index=False)
    combined.to_csv(outdir / "all_peak_fits_with_gates.csv", index=False)
    selection.loc[selection["selected"]].to_csv(
        outdir / "selected_profile_families.csv", index=False
    )

    write_report(
        outdir / "REPORT.md",
        args.scan,
        args.sample_label,
        combined,
        accepted_positions,
        accepted_areas,
        tentative,
    )
    provenance = {
        "scan": int(args.scan),
        "sample_label": args.sample_label,
        "manifest": str(manifest_path),
        "results_root": str(results_root),
        "reporting_rule": (
            "Accepted position = freely resolved strong/detected + delta_BIC>=6 "
            "+ area_SNR>=3 + optimiser success + no centre bound + "
            "position_reportable + q0_reported present; an explicitly "
            "position-only feature may have only its under-resolved FWHM at a "
            "bound. Reference-like components are excluded."
        ),
    }
    (outdir / "report_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Built conservative report in {outdir}")


if __name__ == "__main__":
    main()
