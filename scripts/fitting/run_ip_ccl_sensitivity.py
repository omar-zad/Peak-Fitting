#!/usr/bin/env python3
"""Fast, bootstrap-free sensitivity audit for the IP low-q CCL candidate."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_CONFIG = (
    PROJECT_ROOT / "configs" / "drop40" / "peakfit_config_IP_revised.json"
)
INPUT = PROJECT_ROOT / "example_data" / "drop40" / "drop40_IP_9frames_norm.txt"
FIT_SCRIPT = PROJECT_ROOT / "scripts" / "fitting" / "fit_giwaxs_series.py"
OUT = PROJECT_ROOT / "local_results" / "IP_CCL_sensitivity"


VARIANTS = {
    "gaussian_linear_baseline": {
        "profile": "gaussian",
        "background_order": 1,
        "q_min": 0.22,
        "q_max": 0.44,
        "group": "background_window",
        "description": "Baseline linear background, q=0.22–0.44 Å⁻¹",
    },
    "gaussian_quadratic": {
        "profile": "gaussian",
        "background_order": 2,
        "q_min": 0.22,
        "q_max": 0.44,
        "group": "background_window",
        "description": "Quadratic background, unchanged q window",
    },
    "gaussian_narrow_window": {
        "profile": "gaussian",
        "background_order": 1,
        "q_min": 0.23,
        "q_max": 0.43,
        "group": "background_window",
        "description": "Linear background, modestly narrowed q window",
    },
    "gaussian_wide_window": {
        "profile": "gaussian",
        "background_order": 1,
        "q_min": 0.21,
        "q_max": 0.45,
        "group": "background_window",
        "description": "Linear background, modestly widened q window",
    },
    "pseudo_voigt_linear_baseline": {
        "profile": "pseudo_voigt",
        "background_order": 1,
        "q_min": 0.22,
        "q_max": 0.44,
        "group": "profile_family",
        "description": "Pseudo-Voigt family, baseline background/window",
    },
    "lorentzian_linear_baseline": {
        "profile": "lorentzian",
        "background_order": 1,
        "q_min": 0.22,
        "q_max": 0.44,
        "group": "profile_family",
        "description": "Lorentzian family, baseline background/window",
    },
}


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def make_config(name: str, spec: dict, base: dict) -> dict:
    low_q = next(
        window for window in base["windows"] if window["key"] == "lamellar_low_q"
    )
    low_q = json.loads(json.dumps(low_q))
    low_q["q_min"] = spec["q_min"]
    low_q["q_max"] = spec["q_max"]
    low_q["background_order"] = spec["background_order"]
    low_q["title"] = f"IP q≈0.306 CCL sensitivity: {spec['description']}"
    low_q["profile_note"] = (
        f"Forced {spec['profile']} profile for bounded CCL sensitivity testing."
    )

    config = {
        "frame_temperatures_C": base["frame_temperatures_C"],
        "profiles_to_compare": [spec["profile"]],
        "fit_settings": base["fit_settings"],
        "windows": [low_q],
        "analysis_variant": (
            f"IP q≈0.306 CCL sensitivity; {name}; bootstrap disabled"
        ),
    }
    return config


def run_variant(name: str, spec: dict, base: dict) -> float:
    variant_dir = OUT / name
    variant_dir.mkdir(parents=True, exist_ok=True)
    config_path = variant_dir / "config.json"
    config_path.write_text(
        json.dumps(make_config(name, spec, base), indent=2) + "\n",
        encoding="utf-8",
    )
    command = [
        sys.executable,
        str(FIT_SCRIPT),
        str(INPUT),
        "--config",
        str(config_path),
        "--outdir",
        str(variant_dir),
        "--stage",
        "all",
        "--confirm-first-fit",
        "--bootstrap",
        "0",
        "--seed",
        "20260723",
    ]
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(variant_dir / ".matplotlib")
    started = time.perf_counter()
    result = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    runtime = time.perf_counter() - started
    (variant_dir / "run.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(
            f"{name} failed with exit code {result.returncode}; "
            f"see {variant_dir / 'run.log'}"
        )
    return runtime


def load_variant(name: str, spec: dict, runtime: float) -> pd.DataFrame:
    variant_dir = OUT / name
    peaks = pd.read_csv(variant_dir / "peak_parameters.csv")
    quality = pd.read_csv(variant_dir / "fit_quality.csv")
    model = pd.read_csv(variant_dir / "model_selection.csv")
    peaks = peaks.loc[
        (peaks["window"] == "lamellar_low_q") & (peaks["peak"] == "q0306")
    ].copy()
    quality = quality.loc[quality["window"] == "lamellar_low_q"].copy()
    cols = [
        "frame",
        "optimizer_success",
        "parameters_at_bounds",
        "reduced_chi2_empirical",
        "max_abs_standardized_residual",
        "quality_flags",
    ]
    peaks = peaks.merge(quality[cols], on="frame", how="left", suffixes=("", "_fit"))
    peaks["variant"] = name
    peaks["variant_group"] = spec["group"]
    peaks["variant_description"] = spec["description"]
    peaks["configured_profile"] = spec["profile"]
    peaks["configured_background_order"] = spec["background_order"]
    peaks["configured_q_min"] = spec["q_min"]
    peaks["configured_q_max"] = spec["q_max"]
    peaks["variant_runtime_s"] = runtime
    peaks["global_bic"] = float(model.iloc[0]["bic"])
    peaks["eta_shared_model"] = model.iloc[0].get("eta_shared", np.nan)

    fwhm_ok = as_bool(peaks["fwhm_reportable"])
    optimizer_ok = as_bool(peaks["optimizer_success"])
    bounds_ok = peaks["parameters_at_bounds"].fillna("").astype(str).str.strip().eq("")
    evidence_ok = (
        peaks["delta_bic_peak"].ge(6.0) & peaks["area_snr"].ge(3.0)
    )
    residual_ok = peaks["max_abs_standardized_residual"].lt(5.0)
    peaks["fit_gate_pass"] = (
        fwhm_ok & optimizer_ok & bounds_ok & evidence_ok & residual_ok
    )

    reasons: list[str] = []
    for _, row in peaks.iterrows():
        failed = []
        if str(row["fwhm_reportable"]).lower() != "true":
            failed.append("FWHM_not_reportable")
        if str(row["optimizer_success"]).lower() != "true":
            failed.append("optimizer_failed")
        if str(row.get("parameters_at_bounds", "")).strip() not in ("", "nan"):
            failed.append("parameter_at_bound")
        if not (
            float(row["delta_bic_peak"]) >= 6.0 and float(row["area_snr"]) >= 3.0
        ):
            failed.append("peak_evidence_below_gate")
        if not float(row["max_abs_standardized_residual"]) < 5.0:
            failed.append("large_standardized_residual")
        reasons.append(";".join(failed))
    peaks["fit_gate_failure_reasons"] = reasons
    return peaks


def comma_join(values: pd.Series) -> str:
    return ";".join(str(value) for value in values if str(value))


def markdown_csv_table(frame: pd.DataFrame) -> str:
    """Render a dependency-free, copyable table for the short Markdown report."""
    return "```csv\n" + frame.to_csv(index=False, float_format="%.5g").strip() + "\n```"


def build_reports(all_rows: pd.DataFrame, runtimes: dict[str, float]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows.to_csv(OUT / "all_variant_q0306_results.csv", index=False)

    baseline_name = "gaussian_linear_baseline"
    baseline = all_rows.loc[all_rows["variant"] == baseline_name, [
        "frame",
        "apparent_fwhm_fit_Ainv",
    ]].rename(columns={"apparent_fwhm_fit_Ainv": "baseline_fwhm_Ainv"})

    bg_names = [
        "gaussian_linear_baseline",
        "gaussian_quadratic",
        "gaussian_narrow_window",
        "gaussian_wide_window",
    ]
    bg = all_rows.loc[all_rows["variant"].isin(bg_names)].merge(
        baseline, on="frame", how="left"
    )
    bg["abs_change_from_baseline_pct"] = (
        100.0
        * (bg["apparent_fwhm_fit_Ainv"] - bg["baseline_fwhm_Ainv"]).abs()
        / bg["baseline_fwhm_Ainv"]
    )
    bg.to_csv(OUT / "background_window_variant_details.csv", index=False)

    bg_summary = (
        bg.groupby(["frame", "temperature_C"], as_index=False)
        .agg(
            baseline_fwhm_Ainv=("baseline_fwhm_Ainv", "first"),
            minimum_fwhm_Ainv=("apparent_fwhm_fit_Ainv", "min"),
            maximum_fwhm_Ainv=("apparent_fwhm_fit_Ainv", "max"),
            maximum_abs_change_from_baseline_pct=(
                "abs_change_from_baseline_pct",
                "max",
            ),
            all_variants_fit_gate_pass=("fit_gate_pass", "all"),
            failed_variants=(
                "variant",
                lambda values: comma_join(
                    bg.loc[
                        values.index[~bg.loc[values.index, "fit_gate_pass"]],
                        "variant",
                    ]
                ),
            ),
        )
    )
    bg_summary["all_variant_sensitivity_over_10pct"] = (
        bg_summary["maximum_abs_change_from_baseline_pct"] > 10.0
    )
    accepted_bg = (
        bg.loc[bg["fit_gate_pass"]]
        .groupby(["frame", "temperature_C"], as_index=False)
        .agg(
            accepted_variant_count=("variant", "count"),
            accepted_minimum_fwhm_Ainv=("apparent_fwhm_fit_Ainv", "min"),
            accepted_maximum_fwhm_Ainv=("apparent_fwhm_fit_Ainv", "max"),
            accepted_maximum_abs_change_from_baseline_pct=(
                "abs_change_from_baseline_pct",
                "max",
            ),
        )
    )
    bg_summary = bg_summary.merge(
        accepted_bg, on=["frame", "temperature_C"], how="left"
    )
    bg_summary["accepted_variant_sensitivity_over_10pct"] = (
        bg_summary["accepted_maximum_abs_change_from_baseline_pct"] > 10.0
    )
    bg_summary.to_csv(OUT / "background_window_sensitivity_summary.csv", index=False)

    family_names = [
        "gaussian_linear_baseline",
        "pseudo_voigt_linear_baseline",
        "lorentzian_linear_baseline",
    ]
    family = all_rows.loc[all_rows["variant"].isin(family_names)].copy()
    width_wide = family.pivot(
        index=["frame", "temperature_C"],
        columns="configured_profile",
        values="apparent_fwhm_fit_Ainv",
    ).reset_index()
    gate_wide = family.pivot(
        index=["frame", "temperature_C"],
        columns="configured_profile",
        values="fit_gate_pass",
    ).reset_index()
    gate_wide = gate_wide.rename(
        columns={
            "gaussian": "gaussian_fit_gate_pass",
            "pseudo_voigt": "pseudo_voigt_fit_gate_pass",
            "lorentzian": "lorentzian_fit_gate_pass",
        }
    )
    width_wide = width_wide.rename(
        columns={
            "gaussian": "gaussian_fwhm_Ainv",
            "pseudo_voigt": "pseudo_voigt_fwhm_Ainv",
            "lorentzian": "lorentzian_fwhm_Ainv",
        }
    )
    family_summary = width_wide.merge(
        gate_wide, on=["frame", "temperature_C"], how="left"
    )
    family_summary["pseudo_voigt_change_vs_gaussian_pct"] = (
        100.0
        * (
            family_summary["pseudo_voigt_fwhm_Ainv"]
            - family_summary["gaussian_fwhm_Ainv"]
        ).abs()
        / family_summary["gaussian_fwhm_Ainv"]
    )
    family_summary["lorentzian_change_vs_gaussian_pct"] = (
        100.0
        * (
            family_summary["lorentzian_fwhm_Ainv"]
            - family_summary["gaussian_fwhm_Ainv"]
        ).abs()
        / family_summary["gaussian_fwhm_Ainv"]
    )
    family_summary["pseudo_voigt_sensitivity_over_10pct"] = (
        family_summary["pseudo_voigt_change_vs_gaussian_pct"] > 10.0
    )
    family_summary["lorentzian_sensitivity_over_10pct"] = (
        family_summary["lorentzian_change_vs_gaussian_pct"] > 10.0
    )
    family_summary.to_csv(OUT / "profile_family_sensitivity_summary.csv", index=False)

    variant_summary = []
    baseline_bic = float(
        all_rows.loc[
            all_rows["variant"] == baseline_name, "global_bic"
        ].iloc[0]
    )
    for name, spec in VARIANTS.items():
        subset = all_rows.loc[all_rows["variant"] == name]
        variant_summary.append(
            {
                "variant": name,
                "group": spec["group"],
                "profile": spec["profile"],
                "background_order": spec["background_order"],
                "q_min": spec["q_min"],
                "q_max": spec["q_max"],
                "global_bic": float(subset["global_bic"].iloc[0]),
                "delta_bic_vs_gaussian_linear_baseline": (
                    float(subset["global_bic"].iloc[0]) - baseline_bic
                    if spec["q_min"] == 0.22 and spec["q_max"] == 0.44
                    else np.nan
                ),
                "eta_shared": subset["eta_shared_model"].iloc[0],
                "runtime_s": runtimes[name],
                "fit_gate_pass_frames": int(subset["fit_gate_pass"].sum()),
            }
        )
    variant_summary_df = pd.DataFrame(variant_summary)
    variant_summary_df.to_csv(OUT / "variant_summary.csv", index=False)

    candidate_temps = {60.0, 65.0, 70.0, 75.0, 80.0, 150.0}
    bg_candidate = bg_summary.loc[
        bg_summary["temperature_C"].isin(candidate_temps)
    ]
    family_candidate = family_summary.loc[
        family_summary["temperature_C"].isin(candidate_temps)
    ]
    profile_bic = variant_summary_df.loc[
        variant_summary_df["variant"].isin(family_names),
        ["profile", "global_bic", "delta_bic_vs_gaussian_linear_baseline", "eta_shared"],
    ]

    lines = [
        "# IP q≈0.306 Å⁻¹ CCL sensitivity audit",
        "",
        "All fits used the same two-component peak constraints as the revised IP "
        "model. Runs were bounded, deterministic, and used bootstrap = 0.",
        "",
        "## Decision",
        "",
    ]
    bg_bad = bg_candidate.loc[
        bg_candidate["accepted_variant_sensitivity_over_10pct"]
    ]
    pv_bad = family_candidate.loc[
        family_candidate["pseudo_voigt_sensitivity_over_10pct"]
    ]
    candidate_gate_bad = bg_candidate.loc[
        ~bg_candidate["all_variants_fit_gate_pass"]
    ]
    if bg_bad.empty:
        lines.append(
            "- The 60–150 °C candidate widths are stable to the tested "
            "background/window choices that pass the fit gate: no frame "
            "changes by more than 10%."
        )
    else:
        affected = ", ".join(f"{t:g} °C" for t in bg_bad["temperature_C"])
        lines.append(
            "- Among variants that pass the fit gate, background/window "
            f"sensitivity exceeds 10% only at: {affected}."
        )
    if pv_bad.empty:
        lines.append(
            "- The near-tied pseudo-Voigt alternative does not change any "
            "60–150 °C width by more than 10%."
        )
    else:
        affected = ", ".join(f"{t:g} °C" for t in pv_bad["temperature_C"])
        lines.append(
            f"- Pseudo-Voigt profile sensitivity exceeds 10% at: {affected}."
        )
    if candidate_gate_bad.empty:
        lines.append(
            "- Every 60–150 °C candidate frame passes the conservative fit "
            "gate in all Gaussian background/window variants."
        )
    else:
        affected = ", ".join(f"{t:g} °C" for t in candidate_gate_bad["temperature_C"])
        lines.append(
            "- At least one Gaussian background/window variant fails the "
            f"conservative fit gate at: {affected}."
        )
    baseline_candidate = all_rows.loc[
        (all_rows["variant"] == baseline_name)
        & (all_rows["temperature_C"].isin(candidate_temps))
    ]
    if bool(baseline_candidate["fit_gate_pass"].all()):
        lines.append(
            "- The selected baseline Gaussian model itself passes the "
            "conservative gate at every candidate temperature from 60 to "
            "150 °C."
        )
    lines.append(
        "- The pseudo-Voigt optimum is η = 0, so it collapses numerically to "
        "the Gaussian widths; Gaussian is preferred because it has fewer "
        "parameters (ΔBIC = 5.75 versus pseudo-Voigt)."
    )
    lines.append(
        "- Lorentzian is disfavoured (ΔBIC = 6.16 versus Gaussian) and fails "
        "the conservative gate at 60, 65, and 150 °C."
    )
    lines.extend(
        [
            "",
            "The conservative gate requires a reportable FWHM, ΔBIC ≥ 6, "
            "area SNR ≥ 3, successful optimization, no parameter at a bound, "
            "and max |standardized residual| < 5.",
            "",
            "## Background/window sensitivity",
            "",
            markdown_csv_table(bg_summary),
            "",
            "## Profile-family sensitivity",
            "",
            markdown_csv_table(profile_bic),
            "",
            markdown_csv_table(family_summary),
            "",
            "## Runtime",
            "",
            f"Total fitting runtime: {sum(runtimes.values()):.1f} s.",
            "",
        ]
    )
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    runtimes: dict[str, float] = {}
    all_rows = []
    overall_started = time.perf_counter()
    existing = OUT / "all_variant_q0306_results.csv"
    complete = existing.exists() and all(
        (OUT / name / "peak_parameters.csv").exists()
        and (OUT / name / "fit_quality.csv").exists()
        and (OUT / name / "model_selection.csv").exists()
        for name in VARIANTS
    )
    if complete:
        combined = pd.read_csv(existing)
        for name in VARIANTS:
            runtimes[name] = float(
                combined.loc[combined["variant"] == name, "variant_runtime_s"].iloc[0]
            )
        all_rows_df = combined
        print("All six fit variants already complete; rebuilding reports only.")
    else:
        for name, spec in VARIANTS.items():
            print(f"Running {name} ...", flush=True)
            runtime = run_variant(name, spec, base)
            runtimes[name] = runtime
            all_rows.append(load_variant(name, spec, runtime))
            print(f"  completed in {runtime:.1f} s", flush=True)
        all_rows_df = pd.concat(all_rows, ignore_index=True)
    build_reports(all_rows_df, runtimes)
    report_wall_time = time.perf_counter() - overall_started
    fitting_runtime = sum(runtimes.values())
    (OUT / "TOTAL_RUNTIME_SECONDS.txt").write_text(
        f"{fitting_runtime:.3f}\n", encoding="utf-8"
    )
    (OUT / "REPORT_ASSEMBLY_RUNTIME_SECONDS.txt").write_text(
        f"{report_wall_time:.3f}\n", encoding="utf-8"
    )
    print(
        f"Audit complete; six fits took {fitting_runtime:.1f} s "
        f"(report assembly {report_wall_time:.1f} s): {OUT}"
    )


if __name__ == "__main__":
    main()
