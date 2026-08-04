#!/usr/bin/env python3
"""Measure AgBh/LaB6 line widths in the same FR/IP/OOP geometry as the sample.

The sample line cuts were produced from a 500 x 500 pygix reciprocal-space map
and then averaged into 400 radial q bins between 0.01 and 2.5 A^-1.  This script
repeats that exact map/bin/sector operation for the transmission-geometry
calibrants before fitting their visible rings.  It is a resolution audit, not
an automatic CCL correction.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(PROJECT_ROOT / ".cache" / "matplotlib"),
)

import fabio
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyFAI
import pygix
from scipy.optimize import curve_fit


Q_MIN = 0.01
Q_MAX = 2.5
RADIAL_BINS = 400
MAP_BINS = 500
DQ = (Q_MAX - Q_MIN) / RADIAL_BINS

SECTORS = {
    "IP": (0.0, 20.0),
    "OOP": (70.0, 90.0),
    "FR": (0.0, 90.0),
}

AGBH_Q1 = 0.1076
AGBH_TARGETS = {
    "AgBh_001": 1 * AGBH_Q1,
    "AgBh_003": 3 * AGBH_Q1,
    "AgBh_004": 4 * AGBH_Q1,
    "AgBh_006": 6 * AGBH_Q1,
    "AgBh_007": 7 * AGBH_Q1,
}
LAB6_TARGETS = {
    "LaB6_100": 1.5115,
    "LaB6_110": 2.1375,
}


def read_first_frame(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        return np.asarray(handle["/entry/data/data"][0], dtype=float)


def make_transform(poni: Path, mask: Path) -> pygix.Transform:
    geometry = pyFAI.load(str(poni))
    transform = pygix.Transform()
    transform.load(geometry)
    transform.maskfile = str(mask)
    transform.incident_angle = 0.0
    transform.sample_orientation = 3
    return transform


def reciprocal_map(
    transform: pygix.Transform, image: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    intensity, qxy, qz = transform.transform_reciprocal(
        image,
        npt=(MAP_BINS, MAP_BINS),
        ip_range=(-Q_MAX, 0.05),
        op_range=(0.0, Q_MAX),
        polarization_factor=1.0,
        correctSolidAngle=True,
        method="nearest",
        unit="A",
    )
    intensity = np.asarray(intensity, dtype=float)
    qr = -np.asarray(qxy, dtype=float)
    qz = np.asarray(qz, dtype=float)
    x_order = np.argsort(qr)
    z_order = np.argsort(qz)
    return intensity[z_order, :][:, x_order], qr[x_order], qz[z_order]


def sector_profiles(
    transform: pygix.Transform, image: np.ndarray
) -> pd.DataFrame:
    intensity, qr, qz = reciprocal_map(transform, image)
    qr_grid, qz_grid = np.meshgrid(qr, qz)
    radial_q = np.hypot(qr_grid, qz_grid)
    chi = np.degrees(np.arctan2(qz_grid, qr_grid))
    base = (
        np.isfinite(intensity)
        & (intensity > 1e-8)
        & np.isfinite(radial_q)
        & (qr_grid >= 0)
        & (qz_grid >= 0)
        & (radial_q >= Q_MIN)
        & (radial_q <= Q_MAX)
    )
    edges = np.linspace(Q_MIN, Q_MAX, RADIAL_BINS + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    output: dict[str, np.ndarray] = {"q_Ainv": centres}
    for name, (chi_min, chi_max) in SECTORS.items():
        selected = base & (chi >= chi_min) & (chi <= chi_max)
        sums, _ = np.histogram(
            radial_q[selected], bins=edges, weights=intensity[selected]
        )
        counts, _ = np.histogram(radial_q[selected], bins=edges)
        profile = np.full_like(centres, np.nan)
        profile[counts > 0] = sums[counts > 0] / counts[counts > 0]
        output[name] = profile
        output[f"{name}_pixel_count"] = counts
    return pd.DataFrame(output)


def pseudo_voigt_plus_line(
    q: np.ndarray,
    background: float,
    slope: float,
    height: float,
    center: float,
    fwhm: float,
    eta: float,
) -> np.ndarray:
    relative = (q - center) / fwhm
    gaussian = np.exp(-4.0 * np.log(2.0) * relative**2)
    lorentzian = 1.0 / (1.0 + 4.0 * relative**2)
    return (
        background
        + slope * (q - center)
        + height * ((1.0 - eta) * gaussian + eta * lorentzian)
    )


def fit_ring(
    profile: pd.DataFrame,
    sector: str,
    calibrant: str,
    reflection: str,
    expected_q: float,
) -> tuple[dict[str, object], pd.DataFrame]:
    half_window = 0.035 if expected_q < 1.0 else 0.060
    selected = (
        profile["q_Ainv"].between(expected_q - half_window, expected_q + half_window)
        & profile[sector].notna()
    )
    x = profile.loc[selected, "q_Ainv"].to_numpy(float)
    y = profile.loc[selected, sector].to_numpy(float)
    if len(x) < 9:
        raise RuntimeError(
            f"{calibrant} {reflection} {sector}: only {len(x)} usable points"
        )

    edge_count = max(2, len(x) // 4)
    left = float(np.median(y[:edge_count]))
    right = float(np.median(y[-edge_count:]))
    center_guess = float(x[np.argmax(y)])
    baseline_guess = 0.5 * (left + right)
    slope_guess = (right - left) / max(float(x[-1] - x[0]), 1e-9)
    height_guess = max(float(np.max(y) - baseline_guess), np.ptp(y) / 2, 1e-9)
    center_limit = 0.020 if expected_q < 1.0 else 0.035

    parameters, covariance = curve_fit(
        pseudo_voigt_plus_line,
        x,
        y,
        p0=[
            baseline_guess,
            slope_guess,
            height_guess,
            np.clip(
                center_guess,
                expected_q - 0.9 * center_limit,
                expected_q + 0.9 * center_limit,
            ),
            max(2.0 * DQ, 0.012),
            0.5,
        ],
        bounds=(
            [
                -np.inf,
                -np.inf,
                0.0,
                expected_q - center_limit,
                0.35 * DQ,
                0.0,
            ],
            [
                np.inf,
                np.inf,
                np.inf,
                expected_q + center_limit,
                0.100,
                1.0,
            ],
        ),
        maxfev=30000,
    )
    errors = np.sqrt(np.clip(np.diag(covariance), 0, np.inf))
    fitted = pseudo_voigt_plus_line(x, *parameters)
    residual = y - fitted
    ss_res = float(np.sum(residual**2))
    ss_total = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_total if ss_total > 0 else np.nan
    fwhm = float(parameters[4])
    center_ok = bool(abs(parameters[3] - expected_q) <= DQ)
    fit_ok = bool(np.isfinite(r_squared) and r_squared >= 0.95)
    sampled_ok = bool(fwhm / DQ >= 5.0)
    quantitative_candidate = center_ok and fit_ok and sampled_ok
    row: dict[str, object] = {
        "calibrant": calibrant,
        "reflection": reflection,
        "sector": sector,
        "sector_chi_min_deg": SECTORS[sector][0],
        "sector_chi_max_deg": SECTORS[sector][1],
        "q_expected_Ainv": expected_q,
        "q_fit_Ainv": float(parameters[3]),
        "q_offset_Ainv": float(parameters[3] - expected_q),
        "fwhm_observed_Ainv": fwhm,
        "fwhm_stderr_Ainv": float(errors[4]),
        "eta": float(parameters[5]),
        "points_across_fwhm": fwhm / DQ,
        "r_squared_descriptive": r_squared,
        "center_within_one_q_bin": center_ok,
        "fit_r_squared_at_least_0p95": fit_ok,
        "fwhm_at_least_five_q_bins": sampled_ok,
        "width_use": (
            "quantitative_candidate"
            if quantitative_candidate
            else "resolution_check_only"
        ),
    }
    fitted_data = pd.DataFrame(
        {
            "calibrant": calibrant,
            "reflection": reflection,
            "sector": sector,
            "q_Ainv": x,
            "observed": y,
            "fitted": fitted,
            "residual": residual,
        }
    )
    return row, fitted_data


def run(args: argparse.Namespace) -> None:
    output = args.outdir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    mask_array = np.asarray(fabio.open(str(args.mask)).data)
    agbh = read_first_frame(args.agbh)
    lab6 = read_first_frame(args.lab6)
    if not (mask_array.shape == agbh.shape == lab6.shape):
        raise ValueError("Mask, AgBh, and LaB6 detector shapes do not match")

    transform = make_transform(args.poni, args.mask)
    profiles = {
        "AgBh": sector_profiles(transform, agbh),
        "LaB6": sector_profiles(transform, lab6),
    }
    profiles["AgBh"].to_csv(output / "AgBh_sector_profiles.csv", index=False)
    profiles["LaB6"].to_csv(output / "LaB6_sector_profiles.csv", index=False)

    rows: list[dict[str, object]] = []
    curves: list[pd.DataFrame] = []
    failures: list[str] = []
    for calibrant, targets in (
        ("AgBh", AGBH_TARGETS),
        ("LaB6", LAB6_TARGETS),
    ):
        for reflection, expected in targets.items():
            for sector in SECTORS:
                try:
                    row, fitted = fit_ring(
                        profiles[calibrant],
                        sector,
                        calibrant,
                        reflection,
                        expected,
                    )
                except (RuntimeError, ValueError) as exc:
                    failures.append(str(exc))
                    continue
                rows.append(row)
                curves.append(fitted)

    results = pd.DataFrame(rows)
    fit_curves = pd.concat(curves, ignore_index=True)
    results.to_csv(output / "sector_matched_calibrant_widths.csv", index=False)
    fit_curves.to_csv(output / "sector_matched_calibrant_fits.csv", index=False)

    fig, axes = plt.subplots(
        len(AGBH_TARGETS) + len(LAB6_TARGETS),
        len(SECTORS),
        figsize=(13, 18),
        constrained_layout=True,
    )
    target_rows = [
        ("AgBh", reflection, expected)
        for reflection, expected in AGBH_TARGETS.items()
    ] + [
        ("LaB6", reflection, expected)
        for reflection, expected in LAB6_TARGETS.items()
    ]
    for row_index, (calibrant, reflection, expected) in enumerate(target_rows):
        for column_index, sector in enumerate(SECTORS):
            ax = axes[row_index, column_index]
            subset = fit_curves[
                (fit_curves["calibrant"] == calibrant)
                & (fit_curves["reflection"] == reflection)
                & (fit_curves["sector"] == sector)
            ]
            summary = results[
                (results["calibrant"] == calibrant)
                & (results["reflection"] == reflection)
                & (results["sector"] == sector)
            ]
            if subset.empty:
                ax.text(0.5, 0.5, "fit unavailable", ha="center", va="center")
                ax.set_axis_off()
                continue
            ax.plot(subset["q_Ainv"], subset["observed"], "o", ms=3, label="data")
            ax.plot(subset["q_Ainv"], subset["fitted"], lw=1.8, label="fit")
            ax.axvline(expected, color="0.5", ls=":", lw=1)
            record = summary.iloc[0]
            ax.set_title(
                f"{reflection}, {sector}: "
                f"FWHM={record.fwhm_observed_Ainv:.4f} Å⁻¹ "
                f"({record.points_across_fwhm:.1f} bins)"
            )
            ax.set_xlabel("q (Å⁻¹)")
            ax.set_ylabel("Intensity (a.u.)")
            ax.grid(alpha=0.18)
    fig.suptitle(
        "AgBh and LaB₆ fitted in the sample FR/IP/OOP q-bin geometry",
        fontsize=15,
        y=1.006,
    )
    fig.savefig(output / "sector_matched_calibrant_diagnostics.png", dpi=220)
    plt.close(fig)

    summary_lines = [
        "Sector-matched calibrant width audit",
        "====================================",
        "",
        "Integration: 500 x 500 reciprocal-space map; 400 radial bins from "
        "0.01 to 2.5 A^-1.",
        "Sectors: IP 0-20 deg, OOP 70-90 deg, FR 0-90 deg.",
        "",
        "A calibrant width below five q bins is retained only as a resolution "
        "warning and must not be treated as a precise deconvolution kernel.",
        "The LaB6 widths are local to their q positions and must not be "
        "automatically applied to low-q lamellar peaks.",
        "",
        f"Successful fits: {len(results)}",
        f"Unavailable fits: {len(failures)}",
    ]
    if failures:
        summary_lines.extend(["", "Unavailable:"] + [f"- {item}" for item in failures])
    (output / "README.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )


def parser() -> argparse.ArgumentParser:
    output_default = (
        PROJECT_ROOT / "results" / "calibration" / "sector_matched_calibrants"
    )
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--poni",
        type=Path,
        default=Path("external_data/calibration/AgBH PONI.poni"),
    )
    cli.add_argument(
        "--mask",
        type=Path,
        default=Path("external_data/calibration/AgBH Mask .edf"),
    )
    cli.add_argument(
        "--agbh",
        type=Path,
        default=Path("external_data/calibration/pilatus2-586298.hdf5"),
    )
    cli.add_argument(
        "--lab6",
        type=Path,
        default=Path("external_data/calibration/pilatus2-586303.hdf5"),
    )
    cli.add_argument("--outdir", type=Path, default=output_default)
    return cli


if __name__ == "__main__":
    run(parser().parse_args())
