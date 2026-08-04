#!/usr/bin/env python3
"""Summarise the DROP40 q=0.459/0.480 reciprocal sector sensitivity test."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = (
    PROJECT_ROOT
    / "results"
    / "drop40"
    / "analysis_9frame_fast"
    / "cross_sector_0459_0480_sensitivity"
)


CUTS = {
    "OOP": {
        "test_peak": "q0459_oop_test",
        "test_label": r"added $q\approx0.459$ OOP component",
        "test_plain": "added q~0.459 A^-1 OOP component",
        "established_peak": "q0480_oop",
        "established_label": r"established $q\approx0.481$ OOP peak",
        "other_peak": "q0536_oop",
        "other_label": r"established $q\approx0.535$ OOP peak",
        "baseline_dir": "OOP_revised",
        "baseline_window": "heated_low_q_0480_0536",
        "baseline_peak": "q0480",
    },
    "IP": {
        "test_peak": "q0480_ip_test",
        "test_label": r"added $q\approx0.480$ IP component",
        "test_plain": "added q~0.480 A^-1 IP component",
        "established_peak": "q0459_ip",
        "established_label": r"established $q\approx0.459$ IP peak",
        "other_peak": "q0520_ip",
        "other_label": "existing higher-q IP component",
        "baseline_dir": "IP_revised",
        "baseline_window": "heated_low_q_ip",
        "baseline_peak": "q0461_ip",
    },
}


def _parameter_subset(path: Path, peak: str) -> pd.DataFrame:
    table = pd.read_csv(path)
    return table.loc[table["peak"].eq(peak)].copy()


def build_summary(root: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    baseline_root = root.parent
    for cut, spec in CUTS.items():
        parameters = pd.read_csv(root / cut / "peak_parameters.csv")
        quality = pd.read_csv(root / cut / "fit_quality.csv")
        test = parameters.loc[parameters["peak"].eq(spec["test_peak"])].copy()
        established = parameters.loc[
            parameters["peak"].eq(spec["established_peak"]),
            ["frame", "q0_fit_Ainv"],
        ].rename(columns={"q0_fit_Ainv": "established_q0_Ainv"})

        baseline = pd.read_csv(
            baseline_root / spec["baseline_dir"] / "peak_parameters.csv"
        )
        baseline = baseline.loc[
            baseline["window"].eq(spec["baseline_window"])
            & baseline["peak"].eq(spec["baseline_peak"])
            & baseline["temperature_C"].isin(test["temperature_C"]),
            ["frame", "q0_fit_Ainv"],
        ].rename(columns={"q0_fit_Ainv": "baseline_established_q0_Ainv"})

        merged = test.merge(established, on="frame", how="left")
        merged = merged.merge(baseline, on="frame", how="left")
        merged = merged.merge(
            quality[
                [
                    "frame",
                    "reduced_chi2_empirical",
                    "nrmse_range",
                    "durbin_watson",
                    "max_abs_standardized_residual",
                    "quality_flags",
                ]
            ].rename(columns={"quality_flags": "fit_quality_flags"}),
            on="frame",
            how="left",
        )
        merged.insert(0, "cut", cut)
        merged.insert(3, "tested_counterpart", spec["test_plain"])
        merged["established_q_shift_after_test_Ainv"] = (
            merged["established_q0_Ainv"]
            - merged["baseline_established_q0_Ainv"]
        )
        rows.append(merged)

    combined = pd.concat(rows, ignore_index=True)
    columns = [
        "cut",
        "frame",
        "temperature_C",
        "tested_counterpart",
        "detection_status",
        "delta_bic_peak",
        "area_snr",
        "area_total",
        "area_total_ci95_low",
        "area_total_ci95_high",
        "q0_fit_Ainv",
        "q0_ci95_low_Ainv",
        "q0_ci95_high_Ainv",
        "apparent_fwhm_fit_Ainv",
        "position_reportable",
        "established_q0_Ainv",
        "baseline_established_q0_Ainv",
        "established_q_shift_after_test_Ainv",
        "reduced_chi2_empirical",
        "nrmse_range",
        "durbin_watson",
        "max_abs_standardized_residual",
        "quality_flags",
        "fit_quality_flags",
    ]
    return combined[columns].sort_values(["cut", "temperature_C"])


def plot_overlays(root: Path, output_path: Path) -> None:
    temperatures = [65, 70, 75, 80]
    fig = plt.figure(figsize=(15.0, 9.4))
    grid = fig.add_gridspec(
        5,
        4,
        height_ratios=[3.1, 1.0, 0.34, 3.1, 1.0],
        hspace=0.12,
        wspace=0.24,
    )
    legend_handles = None
    legend_labels = None

    colours = {
        "raw": "#222222",
        "fit": "#D55E00",
        "background": "#777777",
        "test": "#CC79A7",
        "established": "#0072B2",
        "other": "#009E73",
        "residual": "#3366AA",
    }

    for cut_index, (cut, spec) in enumerate(CUTS.items()):
        curves = pd.read_csv(root / cut / "fitted_curves.csv")
        parameters = pd.read_csv(root / cut / "peak_parameters.csv")
        main_row = 0 if cut_index == 0 else 3
        residual_row = 1 if cut_index == 0 else 4
        for column, temperature in enumerate(temperatures):
            data = curves.loc[curves["temperature_C"].eq(temperature)].copy()
            test = parameters.loc[
                parameters["temperature_C"].eq(temperature)
                & parameters["peak"].eq(spec["test_peak"])
            ].iloc[0]

            ax = fig.add_subplot(grid[main_row, column])
            residual_ax = fig.add_subplot(grid[residual_row, column], sharex=ax)
            ax.plot(
                data["q_Ainv"],
                data["observed_intensity"],
                marker="o",
                ms=3.6,
                lw=0.9,
                color=colours["raw"],
                label="raw data",
                zorder=5,
            )
            ax.plot(
                data["q_Ainv"],
                data["total_fit_intensity"],
                lw=2.0,
                color=colours["fit"],
                label="total fit",
                zorder=6,
            )
            ax.plot(
                data["q_Ainv"],
                data["background_intensity"],
                lw=1.1,
                ls="--",
                color=colours["background"],
                label="linear background",
            )
            for role, peak_key in (
                ("test", spec["test_peak"]),
                ("established", spec["established_peak"]),
                ("other", spec["other_peak"]),
            ):
                component = f"component_{peak_key}"
                ax.plot(
                    data["q_Ainv"],
                    data["background_intensity"] + data[component],
                    lw=1.45,
                    ls=":" if role == "test" else "-.",
                    color=colours[role],
                    alpha=0.95,
                    label={
                        "test": "tested cross-sector component",
                        "established": "established main peak",
                        "other": "other retained component",
                    }[role],
                )

            ax.set_title(f"{cut} — {temperature} °C", fontsize=10.5, pad=5)
            ax.text(
                0.03,
                0.04,
                f"test: {test['detection_status'].replace('_', ' ')}\n"
                f"ΔBIC={test['delta_bic_peak']:.1f}, SNR={test['area_snr']:.1f}",
                transform=ax.transAxes,
                fontsize=8.4,
                va="bottom",
                bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 2.5},
            )
            ax.grid(alpha=0.15, lw=0.6)
            ax.tick_params(labelbottom=False)
            if column == 0:
                ax.set_ylabel("Intensity (a.u.)")

            residual_ax.axhspan(-2, 2, color="#999999", alpha=0.13, lw=0)
            residual_ax.axhline(0, color="#555555", lw=0.7)
            residual_ax.plot(
                data["q_Ainv"],
                data["standardized_residual"],
                color=colours["residual"],
                lw=1.1,
            )
            residual_ax.set_ylim(-3.1, 3.1)
            residual_ax.grid(alpha=0.12, lw=0.5)
            residual_ax.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
            if column == 0:
                residual_ax.set_ylabel("residual\n/ noise")

            if legend_handles is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()

    fig.suptitle(
        "DROP40 reciprocal cross-sector sensitivity check\n"
        "Gaussian profiles, linear background, unsmoothed 1D cuts, bootstrap = 0",
        fontsize=14,
        y=0.995,
    )
    fig.legend(
        legend_handles,
        legend_labels,
        loc="lower center",
        ncol=6,
        frameon=False,
        bbox_to_anchor=(0.5, 0.002),
        fontsize=9,
    )
    fig.subplots_adjust(top=0.91, bottom=0.11, left=0.06, right=0.99)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_evidence(summary: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)
    styles = {"OOP": ("#0072B2", "o"), "IP": ("#D55E00", "s")}
    for cut, data in summary.groupby("cut", sort=False):
        colour, marker = styles[cut]
        label = "OOP: test q≈0.459" if cut == "OOP" else "IP: test q≈0.480"
        axes[0].plot(
            data["temperature_C"],
            data["delta_bic_peak"],
            marker=marker,
            color=colour,
            lw=1.6,
            label=label,
        )
        axes[1].plot(
            data["temperature_C"],
            data["area_snr"],
            marker=marker,
            color=colour,
            lw=1.6,
            label=label,
        )
    axes[0].axhline(6, color="#555555", ls="--", lw=1, label="detected threshold")
    axes[0].axhline(2, color="#999999", ls=":", lw=1, label="tentative threshold")
    axes[1].axhline(3, color="#555555", ls="--", lw=1, label="detected threshold")
    axes[1].axhline(2, color="#999999", ls=":", lw=1, label="tentative threshold")
    axes[0].set_ylabel("Leave-one-peak-out ΔBIC")
    axes[1].set_ylabel("Fitted area / standard error")
    for ax in axes:
        ax.set_xlabel("Temperature (°C)")
        ax.set_xticks([65, 70, 75, 80])
        ax.grid(alpha=0.18, lw=0.6)
        ax.legend(frameon=False, fontsize=8.5)
    fig.suptitle("Evidence for the added reciprocal-sector components", fontsize=13)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame, output_path: Path) -> None:
    oop = summary.loc[summary["cut"].eq("OOP")]
    ip = summary.loc[summary["cut"].eq("IP")]
    oop_established = (
        float(oop["established_q0_Ainv"].min()),
        float(oop["established_q0_Ainv"].max()),
    )
    ip_established = (
        float(ip["established_q0_Ainv"].min()),
        float(ip["established_q0_Ainv"].max()),
    )
    oop_shift = float(oop["established_q_shift_after_test_Ainv"].abs().max())
    ip_shift = float(ip["established_q_shift_after_test_Ainv"].abs().max())
    chi2_min = float(summary["reduced_chi2_empirical"].min())
    chi2_max = float(summary["reduced_chi2_empirical"].max())
    max_residual = float(summary["max_abs_standardized_residual"].max())

    lines = [
        "# DROP40 q≈0.459/q≈0.480 reciprocal cross-sector sensitivity check",
        "",
        "## Outcome",
        "",
        "The check supports **sector-selective, separately resolved maxima**, but it does not justify the phrase ‘airtight crystallographic absence’. The OOP test for an added q≈0.459 Å⁻¹ component was negative at every selected temperature. The reciprocal IP test did not find a repeatable q≈0.480 Å⁻¹ counterpart; however, it did label a weak feature nearer q≈0.489–0.491 Å⁻¹ as tentative at 65 and 75 °C. That weak feature was not detected at 70 or 80 °C and never reached the configured detected threshold.",
        "",
        "## Method",
        "",
        "- Frames 197–200: 65, 70, 75 and 80 °C.",
        "- Same q window (0.44–0.57 Å⁻¹), unsmoothed normalized cuts, Gaussian family and linear background as the accepted production fits.",
        "- OOP model retained the established ≈0.481 and ≈0.535 Å⁻¹ components and added a non-overlapping test component around 0.459 Å⁻¹.",
        "- IP model retained the established ≈0.459 and higher-q components and added a non-overlapping test component around 0.480 Å⁻¹.",
        "- Fast covariance run only (bootstrap = 0); this is a model-sensitivity check, not a final uncertainty run.",
        "",
        "## Numerical evidence",
        "",
        f"- OOP added q≈0.459 component: all four frames `not_detected`; ΔBIC {oop['delta_bic_peak'].min():.1f} to {oop['delta_bic_peak'].max():.1f}; maximum area SNR {oop['area_snr'].max():.2f}.",
        f"- IP added q≈0.480 component: `tentative` at 65 and 75 °C (fitted q={ip.loc[ip['detection_status'].eq('tentative'), 'q0_fit_Ainv'].min():.4f}–{ip.loc[ip['detection_status'].eq('tentative'), 'q0_fit_Ainv'].max():.4f} Å⁻¹), and `not_detected` at 70 and 80 °C. Maximum ΔBIC {ip['delta_bic_peak'].max():.1f}; maximum area SNR {ip['area_snr'].max():.2f}; neither reaches the detected thresholds of ΔBIC≥6 and SNR≥3.",
        f"- The established OOP peak remains at {oop_established[0]:.4f}–{oop_established[1]:.4f} Å⁻¹ and shifts by at most {oop_shift:.6f} Å⁻¹ relative to the accepted two-component fit.",
        f"- The established IP peak remains at {ip_established[0]:.4f}–{ip_established[1]:.4f} Å⁻¹ and shifts by at most {ip_shift:.6f} Å⁻¹ relative to the accepted fit.",
        f"- Fit diagnostics remain acceptable: reduced χ² {chi2_min:.2f}–{chi2_max:.2f}; maximum absolute standardised residual {max_residual:.2f}.",
        "- The q step is ≈0.00623 Å⁻¹, so these cuts can exclude a comparable separately resolved counterpart, not an arbitrarily weak or fully overlapping contribution.",
        "",
        "## Dissertation-safe wording",
        "",
        "> Between 70 and 80 °C, a low-q maximum at q≈0.459–0.460 Å⁻¹ was resolved in the IP cut, whereas a distinct maximum at q≈0.481–0.482 Å⁻¹ was resolved in the OOP cut. Reciprocal-component sensitivity fits did not identify a repeatable detected counterpart in the opposite sector, supporting sector-selective scattering within the resolution and detection limits of the 1D cuts. A weak, non-repeatable IP feature near q≈0.489 Å⁻¹ was only tentative and was therefore not assigned as the OOP counterpart.",
        "",
        "Do not replace ‘sector-selective scattering’ with a crystallographic assignment until the peak indexing/assignment and corresponding 2D azimuthal evidence have been agreed with the supervisor.",
        "",
        "## Files to review",
        "",
        "- `cross_sector_sensitivity_overlay.png`: raw data, total fits, components and residuals.",
        "- `cross_sector_sensitivity_evidence.png`: ΔBIC and area-SNR evidence for the added components.",
        "- `cross_sector_sensitivity_summary.csv`: frame-by-frame numerical audit table.",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    root = args.root.resolve()
    summary = build_summary(root)
    summary.to_csv(root / "cross_sector_sensitivity_summary.csv", index=False)
    plot_overlays(root, root / "cross_sector_sensitivity_overlay.png")
    plot_evidence(summary, root / "cross_sector_sensitivity_evidence.png")
    write_report(summary, root / "CROSS_SECTOR_SENSITIVITY_REPORT.md")
    print(f"Wrote summary package to {root}")


if __name__ == "__main__":
    main()
