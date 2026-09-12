#!/usr/bin/env python3
"""Build the explicitly provisional IP q0306 bootstrap/CCL summary."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT = (
    PROJECT_ROOT
    / "results"
    / "drop40"
    / "analysis_9frame_final"
    / "IP_q0306_bootstrap100"
)
PEAKS = OUT / "peak_parameters.csv"
QUALITY = OUT / "fit_quality.csv"

BETA_INST_AINV = 0.012050681085139646
INSTRUMENT_POINTS_ACROSS_FWHM = 1.94
WALL_RUNTIME_SECONDS = 6.58


def corrected_beta(beta: pd.Series) -> pd.Series:
    return np.sqrt(np.maximum(beta**2 - BETA_INST_AINV**2, 0.0))


def main() -> None:
    peaks = pd.read_csv(PEAKS)
    quality = pd.read_csv(QUALITY)
    broad = peaks.loc[peaks["peak"] == "q0306"].copy()
    quality_columns = [
        "frame",
        "reduced_chi2_empirical",
        "max_abs_standardized_residual",
        "optimizer_success",
        "parameters_at_bounds",
        "quality_flags",
        "bootstrap_block_length",
    ]
    broad = broad.merge(
        quality[quality_columns],
        on="frame",
        how="left",
        suffixes=("", "_fit"),
    )

    result = pd.DataFrame(
        {
            "frame": broad["frame"],
            "temperature_C": broad["temperature_C"],
            "q0_Ainv": broad["q0_fit_Ainv"],
            "q0_ci95_low_Ainv": broad["q0_ci95_low_Ainv"],
            "q0_ci95_high_Ainv": broad["q0_ci95_high_Ainv"],
            "observed_fwhm_Ainv": broad["apparent_fwhm_fit_Ainv"],
            "observed_fwhm_ci95_low_Ainv": broad[
                "apparent_fwhm_ci95_low_Ainv"
            ],
            "observed_fwhm_ci95_high_Ainv": broad[
                "apparent_fwhm_ci95_high_Ainv"
            ],
            "bootstrap_successes": broad["bootstrap_successes"],
            "bootstrap_requested": broad["bootstrap_requested"],
            "bootstrap_block_length": broad["bootstrap_block_length"],
            "detection_status": broad["detection_status"],
            "delta_bic_peak": broad["delta_bic_peak"],
            "area_snr": broad["area_snr"],
            "points_across_sample_fwhm": broad["points_across_fwhm"],
            "reduced_chi2_empirical": broad["reduced_chi2_empirical"],
            "max_abs_standardized_residual": broad[
                "max_abs_standardized_residual"
            ],
            "optimizer_success": broad["optimizer_success"],
            "peak_quality_flags": broad["quality_flags"],
            "fit_quality_flags": broad["quality_flags_fit"]
            if "quality_flags_fit" in broad.columns
            else "",
        }
    )

    result["instrument_fwhm_Ainv"] = BETA_INST_AINV
    result["instrument_points_across_fwhm"] = INSTRUMENT_POINTS_ACROSS_FWHM
    result["corrected_fwhm_Ainv"] = corrected_beta(result["observed_fwhm_Ainv"])
    result["corrected_fwhm_ci95_low_Ainv"] = corrected_beta(
        result["observed_fwhm_ci95_low_Ainv"]
    )
    result["corrected_fwhm_ci95_high_Ainv"] = corrected_beta(
        result["observed_fwhm_ci95_high_Ainv"]
    )

    two_pi = 2.0 * np.pi
    result["xi_uncorrected_A"] = two_pi / result["observed_fwhm_Ainv"]
    result["xi_uncorrected_ci95_low_A"] = (
        two_pi / result["observed_fwhm_ci95_high_Ainv"]
    )
    result["xi_uncorrected_ci95_high_A"] = (
        two_pi / result["observed_fwhm_ci95_low_Ainv"]
    )
    result["xi_corrected_provisional_A"] = (
        two_pi / result["corrected_fwhm_Ainv"]
    )
    result["xi_corrected_provisional_ci95_low_A"] = (
        two_pi / result["corrected_fwhm_ci95_high_Ainv"]
    )
    result["xi_corrected_provisional_ci95_high_A"] = (
        two_pi / result["corrected_fwhm_ci95_low_Ainv"]
    )
    for column in [
        "xi_uncorrected_A",
        "xi_uncorrected_ci95_low_A",
        "xi_uncorrected_ci95_high_A",
        "xi_corrected_provisional_A",
        "xi_corrected_provisional_ci95_low_A",
        "xi_corrected_provisional_ci95_high_A",
    ]:
        result[column.replace("_A", "_nm")] = result[column] / 10.0
    result["instrument_correction_effect_on_xi_pct"] = 100.0 * (
        result["xi_corrected_provisional_A"] / result["xi_uncorrected_A"] - 1.0
    )

    result["instrument_correction_status"] = (
        "provisional_Gaussian_quadrature_only"
    )
    result["xi_status"] = (
        "provisional_not_final_dissertation_CCL"
    )
    result["caution"] = (
        "AgBh_003 instrument FWHM spans only 1.94 bins; "
        "resolution correction is approximate"
    )
    result.loc[result["temperature_C"] == 70, "caution"] += (
        "; 70 C also has borderline 10.42% background/window sensitivity"
    )
    result.to_csv(OUT / "IP_q0306_provisional_CCL_table.csv", index=False)

    interval_columns = [
        "frame",
        "temperature_C",
        "q0_Ainv",
        "q0_ci95_low_Ainv",
        "q0_ci95_high_Ainv",
        "observed_fwhm_Ainv",
        "observed_fwhm_ci95_low_Ainv",
        "observed_fwhm_ci95_high_Ainv",
        "bootstrap_successes",
        "bootstrap_requested",
        "detection_status",
        "peak_quality_flags",
        "fit_quality_flags",
    ]
    result[interval_columns].to_csv(
        OUT / "IP_q0306_bootstrap_intervals.csv", index=False
    )

    compact_columns = [
        "temperature_C",
        "q0_Ainv",
        "q0_ci95_low_Ainv",
        "q0_ci95_high_Ainv",
        "observed_fwhm_Ainv",
        "observed_fwhm_ci95_low_Ainv",
        "observed_fwhm_ci95_high_Ainv",
        "xi_uncorrected_nm",
        "xi_corrected_provisional_nm",
        "bootstrap_successes",
        "bootstrap_requested",
    ]
    compact = result[compact_columns].copy()
    compact.columns = [
        "T_C",
        "q0_Ainv",
        "q0_95low",
        "q0_95high",
        "FWHM_Ainv",
        "FWHM_95low",
        "FWHM_95high",
        "xi_uncorrected_nm",
        "xi_corrected_provisional_nm",
        "bootstrap_successes",
        "bootstrap_requested",
    ]

    shoulder = peaks.loc[peaks["peak"] == "q0334"]
    all_success = bool(
        (peaks["bootstrap_successes"] == peaks["bootstrap_requested"]).all()
    )
    q0306_flags = broad["quality_flags"].dropna().astype(str).tolist()
    quality_flags = quality["quality_flags"].dropna().astype(str).tolist()
    lines = [
        "# Final bounded IP q≈0.306 bootstrap result",
        "",
        f"- Wall runtime: {WALL_RUNTIME_SECONDS:.2f} s.",
        "- Scope: IP only; q=0.22–0.44 Å⁻¹; 60, 65, 70, 75, 80, "
        "and 150 °C; two constrained Gaussian components; linear background.",
        "- Bootstrap: 100 residual-bootstrap repeats per frame with fixed "
        "seed 20260723.",
        f"- All peak/frame bootstrap fits successful: {'yes' if all_success else 'no'}.",
        "- q0306 bootstrap successes: "
        + ", ".join(
            f"{temperature:g} °C {int(success)}/100"
            for temperature, success in zip(
                broad["temperature_C"], broad["bootstrap_successes"]
            )
        )
        + ".",
        "- q0306 quality flags: "
        + ("none." if not q0306_flags else "; ".join(q0306_flags) + "."),
        "- Window/frame fit-quality flags: "
        + ("none." if not quality_flags else "; ".join(quality_flags) + "."),
        f"- q0334 shoulder FWHM is under-sampled in all {len(shoulder)} "
        "frames and is not a CCL result; it remains in the model to represent "
        "the shoulder.",
        "",
        "## Provisional q, width, and coherence-length table",
        "",
        "```csv",
        compact.to_csv(index=False, float_format="%.6g").strip(),
        "```",
        "",
        "The uncorrected value is ξ = 2π/βobs. The provisional corrected "
        "value uses Gaussian quadrature, βcorr = √(βobs² − βinst²), with "
        f"βinst = {BETA_INST_AINV:.15f} Å⁻¹ from IP AgBh_003.",
        "For these broad sample peaks, that provisional correction increases "
        "ξ by only "
        f"{result['instrument_correction_effect_on_xi_pct'].min():.2f}–"
        f"{result['instrument_correction_effect_on_xi_pct'].max():.2f}%.",
        "",
        "**Do not present these corrected ξ values as final dissertation CCLs.** "
        "The AgBh width spans only 1.94 q bins, so the instrumental correction "
        "is approximate. In addition, the 70 °C width showed borderline "
        "10.42% model sensitivity. The bootstrap intervals quantify random "
        "residual uncertainty conditional on this selected model; they do not "
        "include those systematic uncertainties.",
        "",
    ]
    (OUT / "FINAL_BOOTSTRAP_REPORT.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    (OUT / "wall_runtime_seconds.txt").write_text(
        f"{WALL_RUNTIME_SECONDS:.2f}\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
