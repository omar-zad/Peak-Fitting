#!/usr/bin/env python3
"""
Reproducible multi-temperature GIWAXS peak fitting.

The workflow compares Gaussian, Lorentzian, and area-normalized pseudo-Voigt
profiles on identical q-ranges and backgrounds.  It then fits the selected
profile to every frame, tests whether each peak is actually supported by the
data, and creates fit/residual panels plus parameter-trend plots.

Example
-------
python scripts/fitting/fit_giwaxs_series.py \
    data/drop40/processed/drop40_FR_9frames_norm.txt \
    --config configs/production/peakfit_config_FR_revised.json \
    --outdir results/drop40/example_run \
    --bootstrap 100

Notes
-----
* The input intensities are fitted exactly as supplied: no smoothing and no
  vertical offsets are applied.
* No pointwise measurement uncertainties were supplied.  Residuals are scaled
  by a robust, frame-local empirical noise estimate solely to make frames
  comparable.  Reported confidence intervals are conditional on the chosen
  peak/background model and do not include systematic uncertainty.
* FWHM values are labelled "reportable" only when the fitted peak spans enough
  sampled q-points and its uncertainty is finite.  Instrumental broadening is
  not removed.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import shutil
import sys
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "base" / "peakfit_config.json"

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(PROJECT_ROOT / ".cache" / "matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.lines import Line2D
from scipy.optimize import least_squares


PROFILE_LABELS = {
    "gaussian": "Gaussian",
    "lorentzian": "Lorentzian",
    "pseudo_voigt": "pseudo-Voigt",
}
PROFILE_COMPLEXITY = {
    "gaussian": 0,
    "lorentzian": 0,
    "pseudo_voigt": 1,
}
COMPONENT_COLORS = [
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
]


def integrate_trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    """Integrate with the NumPy 1.x/2.x compatible trapezoidal API."""
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:
        trapezoid = np.trapz
    return float(trapezoid(y, x))


@dataclass
class CurveFit:
    success: bool
    message: str
    profile: str
    eta_fixed: float | None
    names: list[str]
    x: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    values: dict[str, float]
    stderr: dict[str, float]
    covariance: np.ndarray | None
    condition_number: float
    at_bounds: list[str]
    q: np.ndarray
    observed: np.ndarray
    fitted: np.ndarray
    background: np.ndarray
    components: dict[str, np.ndarray]
    residual: np.ndarray
    standardized_residual: np.ndarray
    noise_sigma: float
    n_parameters: int
    metrics: dict[str, float]


def canonical_frame(value: Any) -> str:
    """Convert 193, '193', or 'Frame 193' to the same frame key."""
    text = str(value).strip()
    matches = re.findall(r"\d+", text)
    if not matches:
        raise ValueError(f"Cannot extract a frame number from column {value!r}")
    return str(int(matches[-1]))


def parse_frame_list(value: str) -> list[str]:
    """Parse a comma-separated, zero-based frame list in acquisition order."""
    frames = [
        canonical_frame(item)
        for item in value.split(",")
        if item.strip()
    ]
    if not frames:
        raise argparse.ArgumentTypeError(
            "--frames requires at least one frame number"
        )
    if len(frames) != len(set(frames)):
        raise argparse.ArgumentTypeError(
            "--frames contains a duplicate frame number"
        )
    return frames


def load_config(path: Path) -> dict[str, Any]:
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {"frame_temperatures_C", "profiles_to_compare", "fit_settings", "windows"}
    missing = sorted(required - set(config))
    if missing:
        raise ValueError(f"Config is missing required keys: {', '.join(missing)}")
    for window in config["windows"]:
        profiles = window.get(
            "profiles_to_compare", config["profiles_to_compare"]
        )
        if not profiles:
            raise ValueError(
                f"{window['key']}: profiles_to_compare cannot be empty"
            )
        unsupported = sorted(set(profiles) - set(PROFILE_LABELS))
        if unsupported:
            raise ValueError(
                f"{window['key']}: unsupported profiles: "
                + ", ".join(unsupported)
            )
        if "eta_fixed" in window:
            eta_fixed = float(window["eta_fixed"])
            if profiles != ["pseudo_voigt"]:
                raise ValueError(
                    f"{window['key']}: eta_fixed requires the window-level "
                    "profiles_to_compare to be ['pseudo_voigt']"
                )
            if not 0.0 <= eta_fixed <= 1.0:
                raise ValueError(
                    f"{window['key']}: eta_fixed must lie between 0 and 1"
                )
        for peak in window["peaks"]:
            if not (peak["q_min"] < peak["q_guess"] < peak["q_max"]):
                raise ValueError(
                    f"{peak['key']}: q_guess must lie strictly inside q_min/q_max"
                )
            if not (
                0 < peak["fwhm_min"]
                <= peak["fwhm_guess"]
                <= peak["fwhm_max"]
            ):
                raise ValueError(
                    f"{peak['key']}: invalid FWHM minimum/guess/maximum"
                )
    return config


def load_data(
    path: Path, temperature_map: dict[str, Any]
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, float], pd.DataFrame]:
    suffix = path.suffix.lower()
    if suffix in {".txt", ".tsv", ".csv", ".dat"}:
        if suffix == ".csv":
            raw = pd.read_csv(path)
        else:
            raw = pd.read_csv(path, sep=None, engine="python")
    elif suffix in {".xlsx", ".xls"}:
        try:
            raw = pd.read_excel(path)
        except ImportError as exc:
            raise RuntimeError(
                "Reading Excel directly requires an Excel engine. Use the paired "
                "tab-delimited .txt file (recommended for reproducibility), or "
                "install openpyxl in the Python environment."
            ) from exc
    else:
        raise ValueError(f"Unsupported input format: {suffix}")

    if raw.shape[1] < 2:
        raise ValueError("Expected q in the first column and at least one frame column")

    q = pd.to_numeric(raw.iloc[:, 0], errors="coerce")
    valid_q = q.notna()
    q_values = q.loc[valid_q].to_numpy(float)
    order = np.argsort(q_values)
    q_values = q_values[order]
    if len(np.unique(q_values)) != len(q_values):
        raise ValueError("q values must be unique")

    frames: dict[str, np.ndarray] = {}
    temperatures: dict[str, float] = {}
    for column in raw.columns[1:]:
        frame = canonical_frame(column)
        if frame in frames:
            raise ValueError(f"Duplicate frame number after normalization: {frame}")
        if frame not in temperature_map:
            raise ValueError(
                f"Frame {frame} has no temperature in frame_temperatures_C"
            )
        y = pd.to_numeric(raw.loc[valid_q, column], errors="coerce").to_numpy(float)
        frames[frame] = y[order]
        temperatures[frame] = float(temperature_map[frame])

    expected = set(map(str, temperature_map))
    missing_frames = sorted(expected - set(frames), key=int)
    if missing_frames:
        raise ValueError(
            "Configured frames are absent from the input: " + ", ".join(missing_frames)
        )

    # Acquisition order is the scientifically meaningful order for in-situ
    # experiments, including long holds with repeated temperatures.
    frame_order = sorted(frames, key=int)
    frames = {frame: frames[frame] for frame in frame_order}
    temperatures = {frame: temperatures[frame] for frame in frame_order}
    return q_values, frames, temperatures, raw


def gaussian_unit(q: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    z = (q - center) / fwhm
    return np.sqrt(4.0 * np.log(2.0) / np.pi) / fwhm * np.exp(
        -4.0 * np.log(2.0) * z * z
    )


def lorentzian_unit(q: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    z = (q - center) / fwhm
    return (2.0 / (np.pi * fwhm)) / (1.0 + 4.0 * z * z)


def profile_unit(
    q: np.ndarray,
    center: float,
    fwhm: float,
    profile: str,
    eta: float | None,
) -> np.ndarray:
    if profile == "gaussian":
        return gaussian_unit(q, center, fwhm)
    if profile == "lorentzian":
        return lorentzian_unit(q, center, fwhm)
    if profile == "pseudo_voigt":
        if eta is None:
            raise ValueError("pseudo-Voigt requires eta")
        return (1.0 - eta) * gaussian_unit(q, center, fwhm) + eta * (
            lorentzian_unit(q, center, fwhm)
        )
    raise ValueError(f"Unknown profile: {profile}")


def robust_noise_sigma(y: np.ndarray) -> float:
    """Robust local white-noise scale estimated from second differences."""
    finite = y[np.isfinite(y)]
    if finite.size < 5:
        return max(float(np.nanstd(finite)), 1e-6)
    second = np.diff(finite, n=2)
    median = np.median(second)
    mad = np.median(np.abs(second - median))
    sigma = mad / (0.6744897501960817 * np.sqrt(6.0))
    dynamic = np.nanpercentile(finite, 95) - np.nanpercentile(finite, 5)
    floor = max(1e-8, dynamic * 1e-5)
    if not np.isfinite(sigma) or sigma <= floor:
        first = np.diff(finite)
        first_mad = np.median(np.abs(first - np.median(first)))
        sigma = first_mad / (0.6744897501960817 * np.sqrt(2.0))
    if not np.isfinite(sigma) or sigma <= floor:
        sigma = max(float(np.nanstd(finite)), floor)
    return max(float(sigma), floor)


def information_metrics(
    residual: np.ndarray,
    standardized_residual: np.ndarray,
    observed: np.ndarray,
    n_parameters: int,
) -> dict[str, float]:
    n = int(observed.size)
    chi2 = float(np.sum(standardized_residual**2))
    raw_sse = float(np.sum(residual**2))
    rmse = float(np.sqrt(raw_sse / max(n, 1)))
    span = float(np.nanmax(observed) - np.nanmin(observed))
    nrmse = rmse / span if span > 0 else np.nan
    scaled_variance = max(chi2 / max(n, 1), np.finfo(float).tiny)
    aic = n * np.log(scaled_variance) + 2.0 * n_parameters
    if n > n_parameters + 1:
        aicc = aic + (
            2.0 * n_parameters * (n_parameters + 1)
            / (n - n_parameters - 1)
        )
    else:
        aicc = np.inf
    bic = n * np.log(scaled_variance) + n_parameters * np.log(max(n, 1))
    total = float(np.sum((observed - np.mean(observed)) ** 2))
    r2 = 1.0 - raw_sse / total if total > 0 else np.nan
    if n > 1:
        dw = float(np.sum(np.diff(residual) ** 2) / max(raw_sse, 1e-30))
        if np.std(residual[:-1]) > 0 and np.std(residual[1:]) > 0:
            lag1 = float(np.corrcoef(residual[:-1], residual[1:])[0, 1])
        else:
            lag1 = np.nan
    else:
        dw = np.nan
        lag1 = np.nan
    return {
        "n_points": float(n),
        "n_parameters": float(n_parameters),
        "chi2_empirical": chi2,
        "reduced_chi2_empirical": chi2 / max(n - n_parameters, 1),
        "rmse_intensity": rmse,
        "nrmse_range": nrmse,
        "aic": float(aic),
        "aicc": float(aicc),
        "bic": float(bic),
        "r_squared_descriptive": float(r2),
        "durbin_watson": dw,
        "residual_lag1": lag1,
        "max_abs_standardized_residual": float(
            np.max(np.abs(standardized_residual))
        ),
    }


def _edge_background_start(
    q: np.ndarray, y: np.ndarray, order: int
) -> tuple[np.ndarray, np.ndarray]:
    midpoint = float((q[0] + q[-1]) / 2.0)
    half_range = max(float((q[-1] - q[0]) / 2.0), 1e-12)
    u = (q - midpoint) / half_range
    edge_count = max(order + 2, int(np.ceil(0.18 * len(q))))
    edge_indices = np.unique(
        np.concatenate(
            [np.arange(min(edge_count, len(q))), np.arange(max(0, len(q) - edge_count), len(q))]
        )
    )
    coefficients = np.polynomial.polynomial.polyfit(
        u[edge_indices], y[edge_indices], deg=order
    )
    return coefficients, u


def _parameter_layout(
    q: np.ndarray,
    y: np.ndarray,
    window: dict[str, Any],
    profile: str,
    eta_fixed: float | None,
    active_peak_keys: set[str],
    fixed_shapes: dict[str, dict[str, float]],
    initial: dict[str, float] | None,
) -> tuple[
    list[str],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, dict[str, float]],
    np.ndarray,
]:
    order = int(window.get("background_order", 1))
    background_start, u = _edge_background_start(q, y, order)
    names: list[str] = []
    x0: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    metadata: dict[str, dict[str, float]] = {}
    initial = initial or {}
    dq = float(np.median(np.diff(q))) if len(q) > 1 else 0.0

    for degree in range(order + 1):
        name = f"background.b{degree}"
        names.append(name)
        x0.append(float(initial.get(name, background_start[degree])))
        lower.append(-np.inf)
        upper.append(np.inf)

    baseline_guess = sum(
        x0[degree] * u**degree for degree in range(order + 1)
    )
    for peak in window["peaks"]:
        key = peak["key"]
        if key not in active_peak_keys:
            continue
        shape = fixed_shapes.get(key, {})
        center_guess = float(
            shape.get(
                "center",
                initial.get(f"{key}.center", peak["q_guess"]),
            )
        )
        fwhm_low = max(float(peak["fwhm_min"]), 1.25 * dq)
        fwhm_high = float(peak["fwhm_max"])
        fwhm_guess = float(
            shape.get(
                "fwhm",
                initial.get(f"{key}.fwhm", peak["fwhm_guess"]),
            )
        )
        center_guess = float(
            np.clip(center_guess, peak["q_min"] + 1e-9, peak["q_max"] - 1e-9)
        )
        fwhm_guess = float(
            np.clip(fwhm_guess, fwhm_low + 1e-9, fwhm_high - 1e-9)
        )
        eta_for_guess = (
            float(eta_fixed)
            if profile == "pseudo_voigt" and eta_fixed is not None
            else float(initial.get("profile.eta", 0.5))
        )
        nearest = int(np.argmin(np.abs(q - center_guess)))
        height_guess = max(float(y[nearest] - baseline_guess[nearest]), 1e-8)
        unit_height = float(
            profile_unit(
                np.array([center_guess]),
                center_guess,
                fwhm_guess,
                profile,
                eta_for_guess if profile == "pseudo_voigt" else None,
            )[0]
        )
        area_guess = max(height_guess / max(unit_height, 1e-12), 1e-10)
        area_guess = float(initial.get(f"{key}.area", area_guess))
        names.append(f"{key}.area")
        x0.append(max(area_guess, 1e-12))
        lower.append(0.0)
        upper.append(np.inf)

        if "center" not in shape:
            names.append(f"{key}.center")
            x0.append(center_guess)
            lower.append(float(peak["q_min"]))
            upper.append(float(peak["q_max"]))
        if "fwhm" not in shape:
            names.append(f"{key}.fwhm")
            x0.append(fwhm_guess)
            lower.append(fwhm_low)
            upper.append(fwhm_high)
        metadata[key] = {
            "center": center_guess,
            "fwhm": fwhm_guess,
        }

    if profile == "pseudo_voigt" and eta_fixed is None:
        names.append("profile.eta")
        x0.append(float(np.clip(initial.get("profile.eta", 0.5), 1e-7, 1 - 1e-7)))
        lower.append(0.0)
        upper.append(1.0)

    return (
        names,
        np.asarray(x0, dtype=float),
        np.asarray(lower, dtype=float),
        np.asarray(upper, dtype=float),
        metadata,
        u,
    )


def _decode_and_evaluate(
    x: np.ndarray,
    names: list[str],
    q: np.ndarray,
    u: np.ndarray,
    window: dict[str, Any],
    profile: str,
    eta_fixed: float | None,
    active_peak_keys: set[str],
    fixed_shapes: dict[str, dict[str, float]],
) -> tuple[
    dict[str, float],
    np.ndarray,
    np.ndarray,
    dict[str, np.ndarray],
]:
    values = {name: float(value) for name, value in zip(names, x)}
    order = int(window.get("background_order", 1))
    background = sum(
        values[f"background.b{degree}"] * u**degree
        for degree in range(order + 1)
    )
    eta = eta_fixed
    if profile == "pseudo_voigt" and eta_fixed is None:
        eta = values["profile.eta"]
    if profile == "pseudo_voigt":
        values["profile.eta"] = float(eta)

    components: dict[str, np.ndarray] = {}
    fitted = background.copy()
    for peak in window["peaks"]:
        key = peak["key"]
        if key not in active_peak_keys:
            continue
        shape = fixed_shapes.get(key, {})
        center = float(
            shape["center"] if "center" in shape else values[f"{key}.center"]
        )
        fwhm = float(
            shape["fwhm"] if "fwhm" in shape else values[f"{key}.fwhm"]
        )
        values[f"{key}.center"] = center
        values[f"{key}.fwhm"] = fwhm
        area = values[f"{key}.area"]
        component = area * profile_unit(q, center, fwhm, profile, eta)
        components[key] = component
        fitted = fitted + component
    return values, fitted, background, components


def _jitter_start(
    x0: np.ndarray,
    names: list[str],
    lower: np.ndarray,
    upper: np.ndarray,
    noise_sigma: float,
    rng: np.random.Generator,
) -> np.ndarray:
    trial = x0.copy()
    for index, name in enumerate(names):
        if name.endswith(".area"):
            trial[index] = max(trial[index] * np.exp(rng.normal(0.0, 0.75)), 1e-12)
        elif name.endswith(".center"):
            span = upper[index] - lower[index]
            trial[index] = trial[index] + rng.normal(0.0, 0.18 * span)
        elif name.endswith(".fwhm"):
            trial[index] = trial[index] * np.exp(rng.normal(0.0, 0.45))
        elif name == "profile.eta":
            trial[index] = rng.uniform(0.03, 0.97)
        elif name.startswith("background."):
            trial[index] = trial[index] + rng.normal(0.0, 0.25 * noise_sigma)
    finite_low = np.isfinite(lower)
    finite_high = np.isfinite(upper)
    trial[finite_low] = np.maximum(
        trial[finite_low], lower[finite_low] + 1e-10
    )
    trial[finite_high] = np.minimum(
        trial[finite_high], upper[finite_high] - 1e-10
    )
    return trial


def fit_curve(
    q: np.ndarray,
    y: np.ndarray,
    window: dict[str, Any],
    profile: str,
    *,
    eta_fixed: float | None = None,
    active_peak_keys: Iterable[str] | None = None,
    fixed_shapes: dict[str, dict[str, float]] | None = None,
    initial: dict[str, float] | None = None,
    n_starts: int = 8,
    rng: np.random.Generator | None = None,
    noise_override: float | None = None,
) -> CurveFit:
    valid = np.isfinite(q) & np.isfinite(y)
    q = np.asarray(q[valid], dtype=float)
    y = np.asarray(y[valid], dtype=float)
    if q.size < 4:
        raise ValueError(f"{window['key']}: fewer than four finite data points")
    active = (
        set(active_peak_keys)
        if active_peak_keys is not None
        else {peak["key"] for peak in window["peaks"]}
    )
    fixed_shapes = fixed_shapes or {}
    noise_sigma = (
        float(noise_override)
        if noise_override is not None
        else robust_noise_sigma(y)
    )
    rng = rng or np.random.default_rng(0)

    names, x0, lower, upper, _, u = _parameter_layout(
        q,
        y,
        window,
        profile,
        eta_fixed,
        active,
        fixed_shapes,
        initial,
    )

    def residual_function(parameters: np.ndarray) -> np.ndarray:
        _, fitted, _, _ = _decode_and_evaluate(
            parameters,
            names,
            q,
            u,
            window,
            profile,
            eta_fixed,
            active,
            fixed_shapes,
        )
        return (fitted - y) / noise_sigma

    starts = [x0]
    for _ in range(max(1, n_starts) - 1):
        starts.append(_jitter_start(x0, names, lower, upper, noise_sigma, rng))

    best = None
    best_objective = np.inf
    messages: list[str] = []
    for start in starts:
        try:
            result = least_squares(
                residual_function,
                start,
                bounds=(lower, upper),
                method="trf",
                x_scale="jac",
                max_nfev=8000,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
        except Exception as exc:  # pragma: no cover - defensive diagnostic path
            messages.append(str(exc))
            continue
        objective = float(np.sum(result.fun**2))
        if np.isfinite(objective) and objective < best_objective:
            best = result
            best_objective = objective
        messages.append(result.message)

    if best is None:
        raise RuntimeError(
            f"{window['key']} {profile}: all multistart fits failed: "
            + " | ".join(messages[:3])
        )

    values, fitted, background, components = _decode_and_evaluate(
        best.x,
        names,
        q,
        u,
        window,
        profile,
        eta_fixed,
        active,
        fixed_shapes,
    )
    residual = y - fitted
    standardized = residual / noise_sigma
    n_parameters = len(names)
    dof = max(len(y) - n_parameters, 1)
    covariance = None
    condition_number = np.inf
    stderr = {name: np.nan for name in names}
    try:
        jtj = best.jac.T @ best.jac
        condition_number = float(np.linalg.cond(jtj))
        covariance = np.linalg.pinv(jtj, rcond=1e-12) * (
            float(np.sum(best.fun**2)) / dof
        )
        diagonal = np.diag(covariance)
        for name, variance in zip(names, diagonal):
            stderr[name] = float(np.sqrt(variance)) if variance >= 0 else np.nan
    except (np.linalg.LinAlgError, ValueError):
        covariance = None

    for key, shape in fixed_shapes.items():
        if key in active:
            if "center" in shape:
                stderr[f"{key}.center"] = np.nan
            if "fwhm" in shape:
                stderr[f"{key}.fwhm"] = np.nan
    if profile == "pseudo_voigt" and eta_fixed is not None:
        stderr["profile.eta"] = np.nan

    at_bounds: list[str] = []
    for name, value, low, high in zip(names, best.x, lower, upper):
        if np.isfinite(low) and np.isfinite(high):
            tolerance = max(1e-7, 0.005 * (high - low))
        else:
            tolerance = 1e-7
        if np.isfinite(low) and value - low <= tolerance:
            at_bounds.append(name)
        elif np.isfinite(high) and high - value <= tolerance:
            at_bounds.append(name)

    metrics = information_metrics(
        residual, standardized, y, n_parameters=n_parameters
    )
    return CurveFit(
        success=bool(best.success),
        message=str(best.message),
        profile=profile,
        eta_fixed=eta_fixed,
        names=names,
        x=best.x,
        lower=lower,
        upper=upper,
        values=values,
        stderr=stderr,
        covariance=covariance,
        condition_number=condition_number,
        at_bounds=at_bounds,
        q=q,
        observed=y,
        fitted=fitted,
        background=background,
        components=components,
        residual=residual,
        standardized_residual=standardized,
        noise_sigma=noise_sigma,
        n_parameters=n_parameters,
        metrics=metrics,
    )


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    if cumulative[-1] <= 0:
        return float(np.median(values))
    return float(values[np.searchsorted(cumulative, cumulative[-1] / 2.0)])


def _global_information(
    fits: dict[str, CurveFit], extra_parameters: int = 0
) -> dict[str, float]:
    standardized = np.concatenate(
        [fit.standardized_residual for fit in fits.values()]
    )
    n = standardized.size
    k = sum(fit.n_parameters for fit in fits.values()) + extra_parameters
    chi2 = float(np.sum(standardized**2))
    scaled_variance = max(chi2 / max(n, 1), np.finfo(float).tiny)
    aic = n * np.log(scaled_variance) + 2.0 * k
    aicc = (
        aic + 2.0 * k * (k + 1) / (n - k - 1)
        if n > k + 1
        else np.inf
    )
    bic = n * np.log(scaled_variance) + k * np.log(max(n, 1))
    return {
        "n_points": float(n),
        "n_parameters": float(k),
        "chi2_empirical": chi2,
        "aic": float(aic),
        "aicc": float(aicc),
        "bic": float(bic),
    }


def compare_profiles(
    q: np.ndarray,
    frame_data: dict[str, np.ndarray],
    window: dict[str, Any],
    profiles: list[str],
    settings: dict[str, Any],
    seed: int,
    eta_override: float | None = None,
) -> tuple[
    str,
    float | None,
    dict[str, CurveFit],
    list[dict[str, Any]],
    dict[str, dict[str, CurveFit]],
]:
    n_starts = int(settings["multistarts_model_selection"])
    all_candidate_fits: dict[str, dict[str, CurveFit]] = {}
    candidate_rows: list[dict[str, Any]] = []
    eta_by_profile: dict[str, float | None] = {}

    for profile in profiles:
        if profile not in PROFILE_LABELS:
            raise ValueError(f"Unsupported profile in config: {profile}")
        if profile != "pseudo_voigt":
            fits: dict[str, CurveFit] = {}
            for frame, y in frame_data.items():
                local_seed = seed + zlib.crc32(
                    f"{window['key']}:{frame}:{profile}".encode()
                )
                fits[frame] = fit_curve(
                    q,
                    y,
                    window,
                    profile,
                    n_starts=n_starts,
                    rng=np.random.default_rng(local_seed),
                )
            metrics = _global_information(fits)
            all_candidate_fits[profile] = fits
            eta_by_profile[profile] = None
            candidate_rows.append(
                {
                    "window": window["key"],
                    "profile": profile,
                    "eta_shared": np.nan,
                    **metrics,
                }
            )
            continue

        if eta_override is not None:
            fixed_eta_fits: dict[str, CurveFit] = {}
            for frame, y in frame_data.items():
                local_seed = seed + zlib.crc32(
                    f"{window['key']}:{frame}:pv-locked:{eta_override:.6f}".encode()
                )
                fixed_eta_fits[frame] = fit_curve(
                    q,
                    y,
                    window,
                    profile,
                    eta_fixed=float(eta_override),
                    n_starts=n_starts,
                    rng=np.random.default_rng(local_seed),
                )
            metrics = _global_information(fixed_eta_fits)
            all_candidate_fits[profile] = fixed_eta_fits
            eta_by_profile[profile] = float(eta_override)
            candidate_rows.append(
                {
                    "window": window["key"],
                    "profile": profile,
                    "eta_shared": float(eta_override),
                    **metrics,
                }
            )
            continue

        free_fits: dict[str, CurveFit] = {}
        eta_values: list[float] = []
        eta_weights: list[float] = []
        for frame, y in frame_data.items():
            local_seed = seed + zlib.crc32(
                f"{window['key']}:{frame}:pv-free".encode()
            )
            fit = fit_curve(
                q,
                y,
                window,
                profile,
                eta_fixed=None,
                n_starts=n_starts,
                rng=np.random.default_rng(local_seed),
            )
            free_fits[frame] = fit
            eta_values.append(fit.values["profile.eta"])
            total_area = sum(
                fit.values.get(f"{peak['key']}.area", 0.0)
                for peak in window["peaks"]
            )
            eta_weights.append(max(total_area / fit.noise_sigma, 1e-6))

        median_eta = _weighted_median(
            np.asarray(eta_values), np.asarray(eta_weights)
        )
        step = float(settings.get("pseudo_voigt_eta_grid_step", 0.1))
        eta_candidates = np.unique(
            np.clip(
                np.asarray(
                    [
                        0.0,
                        0.25,
                        0.5,
                        0.75,
                        1.0,
                        median_eta - 2 * step,
                        median_eta - step,
                        median_eta,
                        median_eta + step,
                        median_eta + 2 * step,
                    ]
                ),
                0.0,
                1.0,
            ).round(6)
        )

        best_eta = None
        best_eta_fits = None
        best_chi2 = np.inf
        for eta in eta_candidates:
            fits = {}
            for frame, y in frame_data.items():
                local_seed = seed + zlib.crc32(
                    f"{window['key']}:{frame}:pv:{eta:.6f}".encode()
                )
                fits[frame] = fit_curve(
                    q,
                    y,
                    window,
                    profile,
                    eta_fixed=float(eta),
                    initial=free_fits[frame].values,
                    n_starts=max(3, n_starts // 2),
                    rng=np.random.default_rng(local_seed),
                    noise_override=free_fits[frame].noise_sigma,
                )
            chi2 = sum(
                fit.metrics["chi2_empirical"] for fit in fits.values()
            )
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_eta = float(eta)
                best_eta_fits = fits

        assert best_eta is not None and best_eta_fits is not None
        metrics = _global_information(best_eta_fits, extra_parameters=1)
        all_candidate_fits[profile] = best_eta_fits
        eta_by_profile[profile] = best_eta
        candidate_rows.append(
            {
                "window": window["key"],
                "profile": profile,
                "eta_shared": best_eta,
                **metrics,
            }
        )

    score = {row["profile"]: row["bic"] for row in candidate_rows}
    selected = min(score, key=score.get)
    if selected == "pseudo_voigt":
        simple = min(
            (profile for profile in score if profile != "pseudo_voigt"),
            key=score.get,
            default=None,
        )
        if simple is not None and score[simple] - score[selected] < 2.0:
            selected = simple

    best_bic = min(score.values())
    for row in candidate_rows:
        row["delta_bic_from_best"] = row["bic"] - best_bic
        row["selected"] = row["profile"] == selected

    return (
        selected,
        eta_by_profile[selected],
        all_candidate_fits[selected],
        candidate_rows,
        all_candidate_fits,
    )


def classify_detection(
    delta_bic: float,
    area_snr: float,
    at_shape_bound: bool,
    settings: dict[str, Any],
) -> str:
    if at_shape_bound:
        if (
            delta_bic >= float(settings["tentative_delta_bic"])
            and area_snr >= float(settings["tentative_area_snr"])
        ):
            return "unresolved"
        return "not_detected"
    if (
        delta_bic >= float(settings["strong_delta_bic"])
        and area_snr >= float(settings["strong_area_snr"])
    ):
        return "strong"
    if (
        delta_bic >= float(settings["detected_delta_bic"])
        and area_snr >= float(settings["detected_area_snr"])
    ):
        return "detected"
    if (
        delta_bic >= float(settings["tentative_delta_bic"])
        and area_snr >= float(settings["tentative_area_snr"])
    ):
        return "tentative"
    return "not_detected"


def test_peak_detection(
    q: np.ndarray,
    y: np.ndarray,
    window: dict[str, Any],
    full_fit: CurveFit,
    profile: str,
    eta_fixed: float | None,
    settings: dict[str, Any],
    seed: int,
) -> dict[str, dict[str, Any]]:
    all_keys = {peak["key"] for peak in window["peaks"]}
    detections: dict[str, dict[str, Any]] = {}
    for peak in window["peaks"]:
        key = peak["key"]
        local_seed = seed + zlib.crc32(
            f"{window['key']}:{key}:detection".encode()
        )
        reduced = fit_curve(
            q,
            y,
            window,
            profile,
            eta_fixed=eta_fixed,
            active_peak_keys=all_keys - {key},
            initial=full_fit.values,
            n_starts=max(4, int(settings["multistarts_final"]) // 2),
            rng=np.random.default_rng(local_seed),
            noise_override=full_fit.noise_sigma,
        )
        delta_bic = reduced.metrics["bic"] - full_fit.metrics["bic"]
        area = full_fit.values[f"{key}.area"]
        area_error = full_fit.stderr.get(f"{key}.area", np.nan)
        area_snr = (
            area / area_error
            if np.isfinite(area_error) and area_error > 0
            else 0.0
        )
        center_at_bound = f"{key}.center" in full_fit.at_bounds
        fwhm_at_bound = f"{key}.fwhm" in full_fit.at_bounds
        # A deliberately under-resolved line may still locate a sampled q
        # maximum even though its width is at the resolution floor.  Explicit
        # position-only candidates therefore use the centre bound, not the
        # FWHM bound, for detection classification.  Their width is suppressed
        # later and must never be interpreted as a line-broadening result.
        at_shape_bound = center_at_bound or (
            fwhm_at_bound and not bool(peak.get("position_only", False))
        )
        status = classify_detection(
            delta_bic, area_snr, at_shape_bound, settings
        )
        detections[key] = {
            "delta_bic_peak": float(delta_bic),
            "area_snr_preliminary": float(area_snr),
            "shape_at_bound_preliminary": bool(at_shape_bound),
            "status": status,
        }
    return detections


def configured_anchor_frame(
    peak: dict[str, Any],
    available_frames: Iterable[str],
    temperatures: dict[str, float],
) -> str:
    """Choose an explicit anchor frame when available, else use temperature."""
    frames = list(available_frames)
    explicit = peak.get("anchor_frame")
    if explicit is not None:
        candidate = canonical_frame(explicit)
        if candidate in frames:
            return candidate
    anchor_temperature = float(peak["anchor_temperature_C"])
    return min(
        frames,
        key=lambda frame: (
            abs(temperatures[frame] - anchor_temperature),
            int(frame),
        ),
    )


def choose_reference_shapes(
    window: dict[str, Any],
    frame_fits: dict[str, CurveFit],
    frame_detections: dict[str, dict[str, dict[str, Any]]],
    temperatures: dict[str, float],
) -> dict[str, dict[str, float]]:
    references: dict[str, dict[str, float]] = {}
    for peak in window["peaks"]:
        key = peak["key"]
        anchor_frame = configured_anchor_frame(
            peak,
            frame_fits,
            temperatures,
        )
        anchor_detection = frame_detections[anchor_frame][key]
        if anchor_detection["status"] not in {"not_detected", "unresolved"}:
            chosen_frame = anchor_frame
        else:
            chosen_frame = max(
                frame_fits,
                key=lambda frame: (
                    frame_detections[frame][key]["delta_bic_peak"],
                    frame_detections[frame][key]["area_snr_preliminary"],
                ),
            )
        fit = frame_fits[chosen_frame]
        references[key] = {
            "center": float(fit.values[f"{key}.center"]),
            "fwhm": float(fit.values[f"{key}.fwhm"]),
            "reference_frame": chosen_frame,
            "reference_temperature_C": float(temperatures[chosen_frame]),
        }
    return references


def refine_reference_shapes_at_anchors(
    q: np.ndarray,
    frame_data: dict[str, np.ndarray],
    window: dict[str, Any],
    preliminary_fits: dict[str, CurveFit],
    preliminary_detections: dict[str, dict[str, dict[str, Any]]],
    references: dict[str, dict[str, float]],
    temperatures: dict[str, float],
    profile: str,
    eta_fixed: float | None,
    settings: dict[str, Any],
    seed: int,
) -> dict[str, dict[str, float]]:
    """Refit each anchor frame while weak peaks from the opposite regime are fixed."""
    anchor_frame_for_peak = {
        peak["key"]: configured_anchor_frame(
            peak,
            frame_data,
            temperatures,
        )
        for peak in window["peaks"]
    }
    anchor_frames = sorted(
        set(anchor_frame_for_peak.values()),
        key=lambda frame: temperatures[frame],
    )
    refined = {key: dict(value) for key, value in references.items()}
    for frame_index, frame in enumerate(anchor_frames):
        fixed_shapes: dict[str, dict[str, float]] = {}
        for peak in window["peaks"]:
            key = peak["key"]
            if anchor_frame_for_peak[key] == frame:
                continue
            if preliminary_detections[frame][key]["status"] in {
                "not_detected",
                "unresolved",
            }:
                fixed_shapes[key] = {
                    "center": refined[key]["center"],
                    "fwhm": refined[key]["fwhm"],
                }
        anchor_fit = fit_curve(
            q,
            frame_data[frame],
            window,
            profile,
            eta_fixed=eta_fixed,
            fixed_shapes=fixed_shapes,
            initial=preliminary_fits[frame].values,
            n_starts=int(settings["multistarts_final"]),
            rng=np.random.default_rng(
                seed
                + frame_index
                + zlib.crc32(f"{window['key']}:{frame}:anchor".encode())
            ),
            noise_override=preliminary_fits[frame].noise_sigma,
        )
        for peak in window["peaks"]:
            key = peak["key"]
            if anchor_frame_for_peak[key] != frame:
                continue
            refined[key] = {
                "center": float(anchor_fit.values[f"{key}.center"]),
                "fwhm": float(anchor_fit.values[f"{key}.fwhm"]),
                "reference_frame": frame,
                "reference_temperature_C": float(temperatures[frame]),
            }
    return refined


def refine_detection_after_stabilization(
    detection: dict[str, Any],
    fit: CurveFit,
    peak: dict[str, Any],
    is_anchor: bool,
    settings: dict[str, Any],
) -> None:
    """Distinguish fixed-shape evidence from a freely resolved line profile."""
    key = peak["key"]
    area = float(fit.values[f"{key}.area"])
    area_error = float(fit.stderr.get(f"{key}.area", np.nan))
    area_snr = (
        area / area_error
        if np.isfinite(area_error) and area_error > 0
        else 0.0
    )
    delta_bic = float(detection["delta_bic_peak"])
    original = detection["status"]
    anchor = bool(is_anchor)
    center_at_bound = f"{key}.center" in fit.at_bounds
    fwhm_at_bound = f"{key}.fwhm" in fit.at_bounds
    final_shape_at_bound = center_at_bound or (
        fwhm_at_bound and not bool(peak.get("position_only", False))
    )
    if (
        anchor
        and not final_shape_at_bound
        and delta_bic >= float(settings["detected_delta_bic"])
        and area_snr >= float(settings["detected_area_snr"])
    ):
        detection["status"] = (
            "strong"
            if delta_bic >= float(settings["strong_delta_bic"])
            and area_snr >= float(settings["strong_area_snr"])
            else "detected"
        )
        return
    if original == "unresolved":
        if (
            delta_bic >= float(settings["detected_delta_bic"])
            and area_snr >= float(settings["detected_area_snr"])
        ):
            detection["status"] = "detected_fixed_shape"
    elif original == "not_detected":
        if (
            anchor
            and delta_bic >= float(settings["strong_delta_bic"])
            and area_snr >= float(settings["strong_area_snr"])
        ):
            detection["status"] = "detected_fixed_shape"
        elif (
            delta_bic >= float(settings["detected_delta_bic"])
            and area_snr >= float(settings["detected_area_snr"])
        ):
            detection["status"] = "tentative_fixed_shape"


def bootstrap_fit(
    fit: CurveFit,
    window: dict[str, Any],
    profile: str,
    eta_fixed: float | None,
    fixed_shapes: dict[str, dict[str, float]],
    count: int,
    seed: int,
) -> tuple[dict[str, list[float]], int, int]:
    if count <= 0:
        return {}, 0, 0
    rng = np.random.default_rng(seed)
    centered_residual = fit.residual - np.mean(fit.residual)
    lag1 = float(fit.metrics.get("residual_lag1", np.nan))
    block_length = 1
    if (
        len(centered_residual) >= 8
        and np.isfinite(lag1)
        and abs(lag1) >= 0.2
    ):
        block_length = max(
            2,
            min(
                int(round(len(centered_residual) ** (1.0 / 3.0))),
                max(2, len(centered_residual) // 3),
            ),
        )
    samples: dict[str, list[float]] = {}
    successes = 0
    for index in range(count):
        if block_length == 1:
            resampled_residual = rng.choice(
                centered_residual,
                size=len(centered_residual),
                replace=True,
            )
        else:
            block_count = int(
                np.ceil(len(centered_residual) / block_length)
            )
            starts = rng.integers(
                0,
                len(centered_residual),
                size=block_count,
            )
            offsets = np.arange(block_length)
            resampled_residual = np.concatenate(
                [
                    centered_residual[
                        (start + offsets) % len(centered_residual)
                    ]
                    for start in starts
                ]
            )[: len(centered_residual)]
        synthetic = fit.fitted + resampled_residual
        try:
            replicate = fit_curve(
                fit.q,
                synthetic,
                window,
                profile,
                eta_fixed=eta_fixed,
                fixed_shapes=fixed_shapes,
                initial=fit.values,
                n_starts=1,
                rng=np.random.default_rng(seed + index + 1),
                noise_override=fit.noise_sigma,
            )
        except Exception:
            continue
        if not replicate.success:
            continue
        successes += 1
        for name, value in replicate.values.items():
            samples.setdefault(name, []).append(float(value))
        for peak in window["peaks"]:
            key = peak["key"]
            component = replicate.components[key]
            samples.setdefault(f"{key}.area_window", []).append(
                integrate_trapezoid(component, replicate.q)
            )
            center = replicate.values[f"{key}.center"]
            samples.setdefault(f"{key}.d_spacing", []).append(
                float(2.0 * np.pi / center)
            )
    return samples, successes, block_length


def confidence_interval(
    name: str,
    value: float,
    stderr: float,
    bootstrap_samples: dict[str, list[float]],
) -> tuple[float, float, str]:
    samples = np.asarray(bootstrap_samples.get(name, []), dtype=float)
    samples = samples[np.isfinite(samples)]
    if samples.size >= 20:
        low, high = np.percentile(samples, [2.5, 97.5])
        return float(low), float(high), "residual_bootstrap_95pct"
    if np.isfinite(stderr) and stderr >= 0:
        return (
            float(value - 1.96 * stderr),
            float(value + 1.96 * stderr),
            "covariance_approx_95pct",
        )
    return np.nan, np.nan, "unavailable"


def fit_quality_flags(fit: CurveFit) -> list[str]:
    flags: list[str] = []
    if not fit.success:
        flags.append("optimizer_not_converged")
    if fit.at_bounds:
        flags.append("parameter_at_bound")
    if not np.isfinite(fit.condition_number) or fit.condition_number > 1e12:
        flags.append("ill_conditioned_covariance")
    lag1 = fit.metrics["residual_lag1"]
    if np.isfinite(lag1) and abs(lag1) > 0.5:
        flags.append("structured_residual")
    if fit.metrics["max_abs_standardized_residual"] > 5.0:
        flags.append("large_standardized_residual")
    return flags


def parameter_rows_for_fit(
    frame: str,
    temperature: float,
    window: dict[str, Any],
    fit: CurveFit,
    detections: dict[str, dict[str, Any]],
    references: dict[str, dict[str, float]],
    bootstrap_samples: dict[str, list[float]],
    bootstrap_successes: int,
    bootstrap_requested: int,
    dq: float,
    settings: dict[str, Any],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    minimum_points = float(
        settings["minimum_points_across_fwhm_for_reporting"]
    )
    fit_flags = fit_quality_flags(fit)
    eta = (
        fit.values.get("profile.eta", np.nan)
        if fit.profile == "pseudo_voigt"
        else (0.0 if fit.profile == "gaussian" else 1.0)
    )
    for peak in window["peaks"]:
        key = peak["key"]
        position_only = bool(peak.get("position_only", False))
        reference_like = bool(peak.get("reference_like_component", False))
        suppress_fwhm = bool(
            peak.get("suppress_fwhm_interpretation", False)
        )
        detection = detections[key]
        status = detection["status"]
        area = float(fit.values[f"{key}.area"])
        center = float(fit.values[f"{key}.center"])
        fwhm = float(fit.values[f"{key}.fwhm"])
        area_error = float(fit.stderr.get(f"{key}.area", np.nan))
        center_error = float(fit.stderr.get(f"{key}.center", np.nan))
        fwhm_error = float(fit.stderr.get(f"{key}.fwhm", np.nan))
        area_window = integrate_trapezoid(fit.components[key], fit.q)
        captured_fraction = (
            area_window / area if area > np.finfo(float).tiny else np.nan
        )
        d_spacing = float(2.0 * np.pi / center)
        d_error = (
            float(2.0 * np.pi * center_error / center**2)
            if np.isfinite(center_error)
            else np.nan
        )
        area_low, area_high, area_ci_method = confidence_interval(
            f"{key}.area", area, area_error, bootstrap_samples
        )
        center_low, center_high, center_ci_method = confidence_interval(
            f"{key}.center", center, center_error, bootstrap_samples
        )
        fwhm_low, fwhm_high, fwhm_ci_method = confidence_interval(
            f"{key}.fwhm", fwhm, fwhm_error, bootstrap_samples
        )
        d_low, d_high, d_ci_method = confidence_interval(
            f"{key}.d_spacing", d_spacing, d_error, bootstrap_samples
        )
        area_window_error = (
            area_error * captured_fraction
            if np.isfinite(area_error) and np.isfinite(captured_fraction)
            else np.nan
        )
        area_window_low, area_window_high, area_window_ci_method = (
            confidence_interval(
                f"{key}.area_window",
                area_window,
                area_window_error,
                bootstrap_samples,
            )
        )
        area_snr = area / area_error if area_error > 0 else np.nan
        points_across = fwhm / dq
        relative_fwhm_error = (
            fwhm_error / fwhm if np.isfinite(fwhm_error) and fwhm > 0 else np.inf
        )
        width_reportable = (
            status in {"strong", "detected", "tentative"}
            and points_across >= minimum_points
            and relative_fwhm_error <= 0.5
            and f"{key}.fwhm" not in fit.at_bounds
            and not suppress_fwhm
            and not reference_like
        )
        position_reportable = (
            status in {"strong", "detected", "tentative"}
            and np.isfinite(center_error)
            and f"{key}.center" not in fit.at_bounds
            and not reference_like
        )
        flags = list(fit_flags)
        if status == "not_detected":
            flags.append("peak_not_detected")
        elif status == "detected_fixed_shape":
            flags.append("peak_detected_with_fixed_shape")
        elif status == "tentative_fixed_shape":
            flags.append("peak_tentative_with_fixed_shape")
        elif status == "unresolved":
            flags.append("peak_unresolved_or_bound")
        elif status == "tentative":
            flags.append("peak_tentative")
        if points_across < minimum_points:
            flags.append("fwhm_undersampled")
        if np.isfinite(captured_fraction) and captured_fraction < 0.9:
            flags.append("profile_tails_outside_window")
        if f"{key}.center" in fit.at_bounds:
            flags.append("peak_shape_at_bound")
        if f"{key}.fwhm" in fit.at_bounds:
            flags.append(
                "fwhm_bound_position_only"
                if position_only
                else "peak_shape_at_bound"
            )
        if position_only:
            flags.append("position_only_width_suppressed")
        if reference_like:
            flags.append("reference_like_component_not_for_assignment")
        if peak.get("provisional_assignment", False):
            flags.append("provisional_unassigned_component")
        if bootstrap_requested and bootstrap_successes < 0.8 * bootstrap_requested:
            flags.append("bootstrap_low_success")

        rows.append(
            {
                "frame": frame,
                "temperature_C": temperature,
                "window": window["key"],
                "window_title": window["title"],
                "peak": key,
                "peak_label": peak["label"],
                "provisional_assignment": bool(
                    peak.get("provisional_assignment", False)
                ),
                "position_only": position_only,
                "reference_like_component": reference_like,
                "suppress_fwhm_interpretation": suppress_fwhm,
                "profile": fit.profile,
                "eta_shared": eta,
                "detection_status": status,
                "delta_bic_peak": detection["delta_bic_peak"],
                "area_snr": area_snr,
                "area_total": area,
                "area_total_ci95_low": area_low,
                "area_total_ci95_high": area_high,
                "area_total_ci_method": area_ci_method,
                "area_in_window": area_window,
                "area_in_window_ci95_low": area_window_low,
                "area_in_window_ci95_high": area_window_high,
                "area_in_window_ci_method": area_window_ci_method,
                "captured_area_fraction": captured_fraction,
                "q0_fit_Ainv": center,
                "q0_ci95_low_Ainv": center_low,
                "q0_ci95_high_Ainv": center_high,
                "q0_ci_method": center_ci_method,
                "q0_reported_Ainv": center if position_reportable else np.nan,
                "d_fit_A": d_spacing,
                "d_ci95_low_A": d_low,
                "d_ci95_high_A": d_high,
                "d_ci_method": d_ci_method,
                "d_reported_A": d_spacing if position_reportable else np.nan,
                "apparent_fwhm_fit_Ainv": fwhm,
                "apparent_fwhm_ci95_low_Ainv": fwhm_low,
                "apparent_fwhm_ci95_high_Ainv": fwhm_high,
                "apparent_fwhm_ci_method": fwhm_ci_method,
                "apparent_fwhm_reported_Ainv": fwhm
                if width_reportable
                else np.nan,
                "points_across_fwhm": points_across,
                "position_reportable": position_reportable,
                "fwhm_reportable": width_reportable,
                "reference_frame_for_weak_fit": references[key][
                    "reference_frame"
                ],
                "reference_temperature_C": references[key][
                    "reference_temperature_C"
                ],
                "bootstrap_requested": bootstrap_requested,
                "bootstrap_successes": bootstrap_successes,
                "quality_flags": ";".join(dict.fromkeys(flags)),
            }
        )
    return rows


def fit_quality_row(
    frame: str,
    temperature: float,
    window: dict[str, Any],
    fit: CurveFit,
) -> dict[str, Any]:
    return {
        "frame": frame,
        "temperature_C": temperature,
        "window": window["key"],
        "window_title": window["title"],
        "profile": fit.profile,
        "eta_shared": fit.values.get("profile.eta", np.nan),
        "noise_sigma_empirical": fit.noise_sigma,
        "condition_number": fit.condition_number,
        "optimizer_success": fit.success,
        "optimizer_message": fit.message,
        "parameters_at_bounds": ";".join(fit.at_bounds),
        "quality_flags": ";".join(fit_quality_flags(fit)),
        **fit.metrics,
    }


def long_rows_for_fit(
    frame: str,
    temperature: float,
    window: dict[str, Any],
    fit: CurveFit,
) -> list[dict[str, Any]]:
    rows = []
    for index, q_value in enumerate(fit.q):
        row: dict[str, Any] = {
            "frame": frame,
            "temperature_C": temperature,
            "window": window["key"],
            "q_Ainv": q_value,
            "observed_intensity": fit.observed[index],
            "total_fit_intensity": fit.fitted[index],
            "background_intensity": fit.background[index],
            "residual_intensity": fit.residual[index],
            "standardized_residual": fit.standardized_residual[index],
        }
        for key, component in fit.components.items():
            row[f"component_{key}"] = component[index]
        rows.append(row)
    return rows


def _plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 9,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def _evenly_spaced_frame_keys(
    frames: list[str], maximum: int
) -> list[str]:
    """Return acquisition-ordered representative keys, including both ends."""
    if len(frames) <= maximum:
        return list(frames)
    indices = np.linspace(0, len(frames) - 1, maximum)
    selected = [frames[int(round(index))] for index in indices]
    return list(dict.fromkeys(selected))


def _plot_window_diagnostic_page(
    window: dict[str, Any],
    frame_fits: dict[str, CurveFit],
    detections: dict[str, dict[str, dict[str, Any]]],
    temperatures: dict[str, float],
    output_path: Path,
    *,
    title_suffix: str = "",
    dpi: int = 240,
) -> None:
    _plot_style()
    frames = list(frame_fits)
    columns = min(4, max(1, len(frames)))
    rows = int(np.ceil(len(frames) / columns))
    fig = plt.figure(
        figsize=(3.8 * columns, 4.35 * rows),
        constrained_layout=False,
    )
    outer = fig.add_gridspec(rows, columns, wspace=0.28, hspace=0.38)
    component_color = {
        peak["key"]: COMPONENT_COLORS[index % len(COMPONENT_COLORS)]
        for index, peak in enumerate(window["peaks"])
    }

    for index, frame in enumerate(frames):
        sub = outer[index // columns, index % columns].subgridspec(
            2, 1, height_ratios=[3.1, 1.0], hspace=0.04
        )
        ax = fig.add_subplot(sub[0])
        residual_ax = fig.add_subplot(sub[1], sharex=ax)
        fit = frame_fits[frame]
        ax.plot(
            fit.q,
            fit.observed,
            "o",
            ms=3.2,
            color="#222222",
            markerfacecolor="white",
            markeredgewidth=0.8,
            label="data" if index == 0 else None,
            zorder=5,
        )
        ax.plot(
            fit.q,
            fit.fitted,
            color="#D55E00",
            lw=1.7,
            label="total fit" if index == 0 else None,
            zorder=6,
        )
        ax.plot(
            fit.q,
            fit.background,
            color="#666666",
            lw=1.1,
            ls="--",
            label="background" if index == 0 else None,
        )
        for peak in window["peaks"]:
            key = peak["key"]
            component_curve = fit.background + fit.components[key]
            status = detections[frame][key]["status"]
            alpha = 1.0 if status != "not_detected" else 0.35
            ax.plot(
                fit.q,
                component_curve,
                color=component_color[key],
                lw=1.15,
                alpha=alpha,
                label=peak["label"] if index == 0 else None,
            )
        eta_text = (
            f", eta={fit.values['profile.eta']:.2f}"
            if fit.profile == "pseudo_voigt"
            else ""
        )
        ax.set_title(
            f"Frame {frame}  |  {temperatures[frame]:g} °C  |  "
            f"{PROFILE_LABELS[fit.profile]}{eta_text}"
        )
        ax.set_ylabel("Intensity (a.u.)")
        ax.grid(alpha=0.16, lw=0.6)
        ax.tick_params(labelbottom=False)

        residual_ax.axhspan(-2, 2, color="#999999", alpha=0.12, lw=0)
        residual_ax.axhline(0, color="#555555", lw=0.7)
        residual_ax.plot(
            fit.q,
            fit.standardized_residual,
            color="#0072B2",
            lw=1.0,
        )
        residual_ax.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
        residual_ax.set_ylabel("r / noise")
        residual_ax.grid(alpha=0.12, lw=0.5)
        robust_limit = max(
            3.0,
            float(np.nanmax(np.abs(fit.standardized_residual))) * 1.08,
        )
        residual_ax.set_ylim(-robust_limit, robust_limit)

    for index in range(len(frames), rows * columns):
        empty_ax = fig.add_subplot(
            outer[index // columns, index % columns]
        )
        empty_ax.axis("off")

    title = (
        f"{window['title']}: observed data, fitted components, and residuals"
    )
    if title_suffix:
        title += f"\n{title_suffix}"
    fig.suptitle(
        title,
        fontsize=13,
        y=0.995,
    )
    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=min(len(labels), 7),
        frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )
    # A single-row diagnostic needs extra headroom so its panel titles do not
    # collide with the figure title. Multi-row pages retain the compact layout.
    fig.subplots_adjust(top=0.87 if rows == 1 else 0.94, bottom=0.09)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_window_diagnostics(
    window: dict[str, Any],
    frame_fits: dict[str, CurveFit],
    detections: dict[str, dict[str, dict[str, Any]]],
    temperatures: dict[str, float],
    output_path: Path,
) -> None:
    """Plot every fit without constructing an unmanageably tall figure."""
    frames = list(frame_fits)
    frames_per_page = 16
    if len(frames) <= frames_per_page:
        _plot_window_diagnostic_page(
            window,
            frame_fits,
            detections,
            temperatures,
            output_path,
        )
        return

    representative_frames = _evenly_spaced_frame_keys(frames, 12)
    _plot_window_diagnostic_page(
        window,
        {frame: frame_fits[frame] for frame in representative_frames},
        detections,
        temperatures,
        output_path,
        title_suffix=(
            f"Representative overview: {len(representative_frames)} of "
            f"{len(frames)} frames; see {output_path.stem}_pages for every frame"
        ),
        dpi=220,
    )

    pages_dir = output_path.parent / f"{output_path.stem}_pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    total_pages = int(np.ceil(len(frames) / frames_per_page))
    for page_index, start in enumerate(
        range(0, len(frames), frames_per_page), start=1
    ):
        page_frames = frames[start : start + frames_per_page]
        _plot_window_diagnostic_page(
            window,
            {frame: frame_fits[frame] for frame in page_frames},
            detections,
            temperatures,
            pages_dir / f"page_{page_index:03d}.png",
            title_suffix=(
                f"All-frame diagnostics, page {page_index} of {total_pages}; "
                f"frames {page_frames[0]}-{page_frames[-1]}"
            ),
            dpi=200,
        )


def _plot_series_overview_page(
    q: np.ndarray,
    frames: dict[str, np.ndarray],
    temperatures: dict[str, float],
    output_path: Path,
    *,
    title_suffix: str = "",
    dpi: int = 240,
) -> None:
    _plot_style()
    frame_count = len(frames)
    columns = min(4, frame_count)
    rows = int(np.ceil(frame_count / columns))
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(3.8 * columns, 3.0 * rows),
        sharex=True,
        constrained_layout=True,
        squeeze=False,
    )
    for ax, (frame, y) in zip(axes.flat, frames.items()):
        valid = np.isfinite(y)
        ax.plot(q[valid], y[valid], color="#0072B2", lw=1.15)
        ax.set_title(f"Frame {frame}  |  {temperatures[frame]:g} °C")
        ax.grid(alpha=0.18, lw=0.6)
        ax.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
        ax.set_ylabel("Intensity (a.u.)")
    for ax in axes.flat[frame_count:]:
        ax.set_visible(False)
    title = "Unsmoothed, un-offset GIWAXS line cuts"
    if title_suffix:
        title += f"\n{title_suffix}"
    fig.suptitle(title, fontsize=13)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_series_overview(
    q: np.ndarray,
    frames: dict[str, np.ndarray],
    temperatures: dict[str, float],
    output_path: Path,
) -> None:
    """Plot a compact overview plus paginated raw curves for long series."""
    frame_keys = list(frames)
    frames_per_page = 20
    if len(frame_keys) <= frames_per_page:
        _plot_series_overview_page(
            q, frames, temperatures, output_path
        )
        return

    representative_frames = _evenly_spaced_frame_keys(frame_keys, 12)
    _plot_series_overview_page(
        q,
        {frame: frames[frame] for frame in representative_frames},
        temperatures,
        output_path,
        title_suffix=(
            f"Representative overview: {len(representative_frames)} of "
            f"{len(frame_keys)} frames; see {output_path.stem}_pages for every frame"
        ),
        dpi=220,
    )

    pages_dir = output_path.parent / f"{output_path.stem}_pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    total_pages = int(np.ceil(len(frame_keys) / frames_per_page))
    for page_index, start in enumerate(
        range(0, len(frame_keys), frames_per_page), start=1
    ):
        page_frames = frame_keys[start : start + frames_per_page]
        _plot_series_overview_page(
            q,
            {frame: frames[frame] for frame in page_frames},
            temperatures,
            pages_dir / f"page_{page_index:03d}.png",
            title_suffix=(
                f"All-frame raw curves, page {page_index} of {total_pages}; "
                f"frames {page_frames[0]}-{page_frames[-1]}"
            ),
            dpi=200,
        )


def _peak_order(config: dict[str, Any]) -> list[tuple[str, str, str]]:
    return [
        (window["key"], peak["key"], peak["label"])
        for window in config["windows"]
        for peak in window["peaks"]
    ]


def plot_parameter_small_multiples(
    parameters: pd.DataFrame,
    config: dict[str, Any],
    output_path: Path,
    metric: str,
    x_axis: str = "temperature",
) -> None:
    _plot_style()
    peak_order = _peak_order(config)
    columns = 6 if len(peak_order) > 16 else 4
    rows = int(np.ceil(len(peak_order) / columns))
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(18.5 if columns == 6 else 14.8, 2.85 * rows),
        constrained_layout=True,
        squeeze=False,
    )
    settings = {
        "area": {
            "value": "area_in_window",
            "low": "area_in_window_ci95_low",
            "high": "area_in_window_ci95_high",
            "ylabel": "Integrated intensity\n(in fitted q-window)",
            "title": "Peak area",
        },
        "position": {
            "value": "q0_reported_Ainv",
            "low": "q0_ci95_low_Ainv",
            "high": "q0_ci95_high_Ainv",
            "ylabel": r"$q_0$ ($\mathrm{\AA}^{-1}$)",
            "title": "Reportable peak position",
        },
        "fwhm": {
            "value": "apparent_fwhm_reported_Ainv",
            "low": "apparent_fwhm_ci95_low_Ainv",
            "high": "apparent_fwhm_ci95_high_Ainv",
            "ylabel": r"Apparent FWHM ($\mathrm{\AA}^{-1}$)",
            "title": "Reportable apparent FWHM",
        },
    }[metric]

    for ax, (window, peak, label) in zip(axes.flat, peak_order):
        data = parameters[
            (parameters["window"] == window) & (parameters["peak"] == peak)
        ].copy()
        data["frame_number"] = data["frame"].map(
            lambda value: int(canonical_frame(value))
        )
        if x_axis == "frame":
            data = data.sort_values("frame_number")
            x_value = data["frame_number"].to_numpy(float)
            x_label = "Acquisition frame (zero-indexed)"
        elif x_axis == "temperature":
            data = data.sort_values("temperature_C")
            x_value = data["temperature_C"].to_numpy(float)
            x_label = "Nominal temperature (°C)"
        else:
            raise ValueError(f"Unsupported parameter-plot x-axis: {x_axis}")
        value = data[settings["value"]].to_numpy(float)
        low = data[settings["low"]].to_numpy(float)
        high = data[settings["high"]].to_numpy(float)
        valid = np.isfinite(value)
        lower_error = np.maximum(value - low, 0)
        upper_error = np.maximum(high - value, 0)
        error = np.vstack([lower_error, upper_error])
        ax.errorbar(
            x_value[valid],
            value[valid],
            yerr=error[:, valid] if valid.any() else None,
            color="#0072B2",
            marker="o",
            ms=4,
            lw=1.1,
            capsize=2,
        )
        if metric == "area":
            tentative = data["detection_status"].isin(
                ["tentative", "tentative_fixed_shape"]
            ).to_numpy()
            fixed_shape = data["detection_status"].eq(
                "detected_fixed_shape"
            ).to_numpy()
            absent = data["detection_status"].eq("not_detected").to_numpy()
            ax.scatter(
                x_value[tentative],
                value[tentative],
                marker="^",
                facecolors="white",
                edgecolors="#E69F00",
                zorder=5,
                s=28,
            )
            ax.scatter(
                x_value[fixed_shape],
                value[fixed_shape],
                marker="s",
                facecolors="white",
                edgecolors="#009E73",
                zorder=5,
                s=26,
            )
            upper_limit = high.copy()
            upper_limit[~np.isfinite(upper_limit)] = value[
                ~np.isfinite(upper_limit)
            ]
            ax.scatter(
                x_value[absent],
                upper_limit[absent],
                marker="v",
                facecolors="white",
                edgecolors="#777777",
                zorder=5,
                s=26,
            )
            ax.set_ylim(bottom=0)
        if not valid.any():
            ax.text(
                0.5,
                0.5,
                "No reportable values",
                transform=ax.transAxes,
                ha="center",
                va="center",
                color="#666666",
                fontsize=8,
            )
        ax.set_title(label)
        ax.grid(alpha=0.18, lw=0.6)
        ax.set_xlabel(x_label)
        ax.set_ylabel(settings["ylabel"])

    for ax in axes.flat[len(peak_order) :]:
        ax.set_visible(False)
    axis_title = (
        "acquisition frame" if x_axis == "frame" else "temperature"
    )
    fig.suptitle(f"{settings['title']} versus {axis_title}", fontsize=13)
    if metric == "area":
        status_handles = [
            Line2D(
                [0],
                [0],
                color="#0072B2",
                marker="o",
                ms=5,
                lw=1.1,
                label="Freely resolved / detected",
            ),
            Line2D(
                [0],
                [0],
                color="none",
                marker="s",
                markerfacecolor="white",
                markeredgecolor="#009E73",
                ms=6,
                label="Detected; shape fixed",
            ),
            Line2D(
                [0],
                [0],
                color="none",
                marker="^",
                markerfacecolor="white",
                markeredgecolor="#E69F00",
                ms=6,
                label="Tentative",
            ),
            Line2D(
                [0],
                [0],
                color="none",
                marker="v",
                markerfacecolor="white",
                markeredgecolor="#777777",
                ms=6,
                label="Not detected; upper 95% CI",
            ),
        ]
        fig.legend(
            handles=status_handles,
            loc="lower center",
            ncol=4,
            frameon=False,
            bbox_to_anchor=(0.5, -0.015),
        )
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_model_selection(
    selection: pd.DataFrame, output_path: Path
) -> None:
    _plot_style()
    windows = list(dict.fromkeys(selection["window"]))
    profiles = [profile for profile in PROFILE_LABELS if profile in set(selection["profile"])]
    x = np.arange(len(windows))
    width = 0.22
    fig, ax = plt.subplots(figsize=(13.5, 5.2), constrained_layout=True)
    for index, profile in enumerate(profiles):
        values = []
        labels = []
        for window in windows:
            matches = selection[
                (selection["window"] == window)
                & (selection["profile"] == profile)
            ]
            if matches.empty:
                values.append(np.nan)
                labels.append(np.nan)
                continue
            row = matches.iloc[0]
            value = float(row["delta_bic_from_best"])
            values.append(min(value, 50.0))
            labels.append(value)
        positions = x + (index - (len(profiles) - 1) / 2.0) * width
        bars = ax.bar(
            positions,
            values,
            width=width,
            label=PROFILE_LABELS[profile],
            color=COMPONENT_COLORS[index],
            alpha=0.82,
        )
        for bar, actual in zip(bars, labels):
            if np.isfinite(actual) and actual > 50:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    49,
                    f"{actual:.0f}+",
                    rotation=90,
                    ha="center",
                    va="top",
                    fontsize=7,
                    color="white",
                )
    ax.axhline(2, color="#555555", lw=0.8, ls="--")
    ax.axhline(6, color="#888888", lw=0.8, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels(windows, rotation=25, ha="right")
    ax.set_ylabel("ΔBIC from best profile (clipped at 50)")
    ax.set_title(
        "Line-shape model comparison across each configured temperature subset"
    )
    ax.legend(frameon=False, ncol=len(profiles))
    ax.grid(axis="y", alpha=0.18)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def _plot_qc_heatmap_page(
    quality: pd.DataFrame,
    output_path: Path,
    *,
    title_suffix: str = "",
    dpi: int = 240,
) -> None:
    _plot_style()
    working = quality.copy()
    working["frame_key"] = working["frame"].map(canonical_frame)
    frame_order = list(dict.fromkeys(working["frame_key"]))
    temperature_by_frame = (
        working.drop_duplicates("frame_key")
        .set_index("frame_key")["temperature_C"]
        .to_dict()
    )
    pivot = working.pivot(
        index="window",
        columns="frame_key",
        values="reduced_chi2_empirical",
    ).reindex(columns=frame_order)
    pivot.columns.name = None
    pivot.index.name = None
    x_labels = [
        f"f{frame}\n{float(temperature_by_frame[frame]):g} °C"
        for frame in frame_order
    ]
    fig, ax = plt.subplots(
        figsize=(max(10.5, 0.65 * len(frame_order)), 5.5),
        constrained_layout=True,
    )
    values = pivot.to_numpy(float)
    finite = values[np.isfinite(values)]
    vmax = float(np.percentile(finite, 90)) if finite.size else 1.0
    image = ax.imshow(
        values,
        aspect="auto",
        cmap="viridis",
        vmin=0,
        vmax=max(vmax, 1.0),
    )
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(x_labels)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("Frame and measured temperature")
    title = (
        "Empirical reduced chi-square "
        "(noise estimated locally; use comparatively)"
    )
    if title_suffix:
        title += f"\n{title_suffix}"
    ax.set_title(title)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            if np.isfinite(value):
                text_color = "white" if value > 0.55 * max(vmax, 1.0) else "black"
                ax.text(
                    column,
                    row,
                    f"{value:.1f}",
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize=8,
                )
    colorbar = fig.colorbar(image, ax=ax, shrink=0.82)
    colorbar.set_label("Empirical reduced χ²")
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_qc_heatmap(quality: pd.DataFrame, output_path: Path) -> None:
    """Plot a readable QC overview and paginate long acquisitions."""
    working = quality.copy()
    working["frame_key"] = working["frame"].map(canonical_frame)
    frame_order = list(dict.fromkeys(working["frame_key"]))
    frames_per_page = 25
    if len(frame_order) <= frames_per_page:
        _plot_qc_heatmap_page(quality, output_path)
        return

    representative_frames = _evenly_spaced_frame_keys(frame_order, 20)
    representative = working[
        working["frame_key"].isin(representative_frames)
    ].drop(columns="frame_key")
    _plot_qc_heatmap_page(
        representative,
        output_path,
        title_suffix=(
            f"Representative overview: {len(representative_frames)} of "
            f"{len(frame_order)} frames; see {output_path.stem}_pages for every frame"
        ),
        dpi=220,
    )

    pages_dir = output_path.parent / f"{output_path.stem}_pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    total_pages = int(np.ceil(len(frame_order) / frames_per_page))
    for page_index, start in enumerate(
        range(0, len(frame_order), frames_per_page), start=1
    ):
        page_frames = frame_order[start : start + frames_per_page]
        page = working[
            working["frame_key"].isin(page_frames)
        ].drop(columns="frame_key")
        _plot_qc_heatmap_page(
            page,
            pages_dir / f"page_{page_index:03d}.png",
            title_suffix=(
                f"All-frame QC, page {page_index} of {total_pages}; "
                f"frames {page_frames[0]}-{page_frames[-1]}"
            ),
            dpi=220,
        )


def plot_first_frame_diagnostic(
    config: dict[str, Any],
    fits_by_window: dict[str, dict[str, CurveFit]],
    parameters: pd.DataFrame,
    quality: pd.DataFrame,
    frame: str,
    temperature: float,
    output_path: Path,
) -> None:
    """Create the single review figure used to unlock the all-frame stage."""
    _plot_style()
    windows = [
        window
        for window in config["windows"]
        if window["key"] in fits_by_window
        and frame in fits_by_window[window["key"]]
    ]
    columns = 2
    rows = int(np.ceil(len(windows) / columns))
    fig = plt.figure(figsize=(19.5, 4.8 * rows))
    outer = fig.add_gridspec(
        rows,
        columns,
        hspace=0.34,
        wspace=0.18,
    )

    for index, window in enumerate(windows):
        row, column = divmod(index, columns)
        inner = outer[row, column].subgridspec(
            2,
            2,
            height_ratios=[3.0, 1.0],
            width_ratios=[3.6, 1.7],
            hspace=0.05,
            wspace=0.08,
        )
        fit = fits_by_window[window["key"]][frame]
        main_ax = fig.add_subplot(inner[0, 0])
        residual_ax = fig.add_subplot(inner[1, 0], sharex=main_ax)
        text_ax = fig.add_subplot(inner[:, 1])
        text_ax.axis("off")

        main_ax.plot(
            fit.q,
            fit.observed,
            color="#111111",
            lw=1.25,
            label="Raw data",
        )
        main_ax.plot(
            fit.q,
            fit.fitted,
            color="#D55E00",
            lw=1.8,
            label="Total fit",
        )
        main_ax.plot(
            fit.q,
            fit.background,
            color="#777777",
            lw=1.0,
            ls="--",
            label="Local background",
        )
        for peak_index, peak in enumerate(window["peaks"]):
            color = COMPONENT_COLORS[peak_index % len(COMPONENT_COLORS)]
            main_ax.plot(
                fit.q,
                fit.background + fit.components[peak["key"]],
                color=color,
                lw=1.0,
                ls=(0, (3, 2)),
                alpha=0.9,
            )
        main_ax.set_title(
            f"{window['title']}  |  {PROFILE_LABELS[fit.profile]}",
            fontsize=10.5,
            loc="left",
        )
        main_ax.set_ylabel("Intensity (a.u.)")
        main_ax.grid(alpha=0.16, lw=0.55)
        main_ax.tick_params(labelbottom=False)

        residual_ax.axhspan(-2, 2, color="#999999", alpha=0.12, lw=0)
        residual_ax.axhline(0, color="#555555", lw=0.7)
        residual_ax.axhline(3, color="#E69F00", lw=0.6, ls=":")
        residual_ax.axhline(-3, color="#E69F00", lw=0.6, ls=":")
        residual_ax.plot(
            fit.q,
            fit.standardized_residual,
            color="#0072B2",
            lw=1.0,
        )
        residual_ax.set_xlabel(r"$q$ ($\mathrm{\AA}^{-1}$)")
        residual_ax.set_ylabel("r/noise")
        residual_ax.grid(alpha=0.12, lw=0.5)
        residual_limit = max(
            3.5,
            1.08
            * float(
                np.nanmax(
                    np.abs(fit.standardized_residual),
                )
            ),
        )
        residual_ax.set_ylim(-residual_limit, residual_limit)

        quality_row = quality[
            (quality["frame"].astype(str) == str(frame))
            & (quality["window"] == window["key"])
        ].iloc[0]
        eta_text = (
            f"; eta={fit.eta_fixed:.3f}"
            if fit.profile == "pseudo_voigt" and fit.eta_fixed is not None
            else ""
        )
        text_ax.text(
            0.0,
            1.0,
            f"Fit quality\n"
            f"profile: {PROFILE_LABELS[fit.profile]}{eta_text}\n"
            f"NRMSE: {quality_row['nrmse_range']:.4f}\n"
            f"empirical reduced chi2: "
            f"{quality_row['reduced_chi2_empirical']:.2f}\n"
            f"R2 (descriptive): "
            f"{quality_row['r_squared_descriptive']:.4f}\n"
            f"Durbin-Watson: {quality_row['durbin_watson']:.2f}\n"
            f"QC: {quality_row['quality_flags'] or 'none'}",
            ha="left",
            va="top",
            fontsize=8.0,
            linespacing=1.25,
        )
        peak_rows = parameters[
            (parameters["frame"].astype(str) == str(frame))
            & (parameters["window"] == window["key"])
        ].set_index("peak")
        y_text = 0.56
        for peak_index, peak in enumerate(window["peaks"]):
            result = peak_rows.loc[peak["key"]]
            color = COMPONENT_COLORS[peak_index % len(COMPONENT_COLORS)]
            q0_text = (
                f"{result['q0_fit_Ainv']:.4f}"
                if bool(result["position_reportable"])
                else f"not reported (fit {result['q0_fit_Ainv']:.4f})"
            )
            fwhm_text = (
                f"{result['apparent_fwhm_fit_Ainv']:.4f}"
                if bool(result["fwhm_reportable"])
                else (
                    f"not reported "
                    f"(fit {result['apparent_fwhm_fit_Ainv']:.4f}; "
                    f"{result['points_across_fwhm']:.2f} steps)"
                )
            )
            text_ax.text(
                0.0,
                y_text,
                f"{peak['label']} [{result['detection_status']}]\n"
                f"q0={q0_text}; Awin={result['area_in_window']:.4g}\n"
                f"FWHM={fwhm_text}",
                ha="left",
                va="top",
                fontsize=7.6,
                color=color,
                linespacing=1.15,
            )
            y_text -= 0.13
        text_ax.text(
            0.0,
            0.01,
            "Dashed coloured curves = individual peak + fitted background.",
            ha="left",
            va="bottom",
            fontsize=7.2,
            color="#666666",
            wrap=True,
        )

        if index == 0:
            handles, labels = main_ax.get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="upper center",
                ncol=3,
                frameon=False,
                bbox_to_anchor=(0.5, 0.963),
            )

    for index in range(len(windows), rows * columns):
        row, column = divmod(index, columns)
        empty_ax = fig.add_subplot(outer[row, column])
        empty_ax.axis("off")

    fig.suptitle(
        f"Stage 1 fit gate: frame {frame}, {temperature:g} °C",
        fontsize=15,
        y=0.995,
    )
    fig.subplots_adjust(top=0.915)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    output_path: Path,
    selection: pd.DataFrame,
    parameters: pd.DataFrame,
    quality: pd.DataFrame,
    dq: float,
    bootstrap: int,
    minimum_points_across_fwhm: float,
) -> None:
    selected = selection[selection["selected"]].copy()
    status_counts = (
        parameters["detection_status"].value_counts().reindex(
            [
                "strong",
                "detected",
                "detected_fixed_shape",
                "tentative",
                "tentative_fixed_shape",
                "unresolved",
                "not_detected",
            ],
            fill_value=0,
        )
    )
    quality_flag_count = int(quality["quality_flags"].fillna("").ne("").sum())
    width_count = int(parameters["fwhm_reportable"].sum())
    lines = [
        "GIWAXS peak-fit run summary",
        "============================",
        "",
        f"q sampling interval: {dq:.7f} A^-1",
        f"Moving-block residual-bootstrap replicates requested per "
        f"window/frame: {bootstrap}",
        "",
        "Selected profile by fitting window:",
    ]
    for row in selected.itertuples(index=False):
        eta = (
            f", shared eta={row.eta_shared:.3f}"
            if row.profile == "pseudo_voigt"
            else ""
        )
        lines.append(
            f"  - {row.window}: {PROFILE_LABELS[row.profile]}{eta}"
        )
    lines.extend(
        [
            "",
            "Peak detection classifications:",
            *[
                f"  - {status}: {int(count)}"
                for status, count in status_counts.items()
            ],
            "",
            f"Reportable apparent FWHM values: {width_count} / {len(parameters)}",
            f"Minimum q steps across FWHM for width reporting: "
            f"{minimum_points_across_fwhm:g}",
            f"Window/frame fits carrying at least one QC flag: "
            f"{quality_flag_count} / {len(quality)}",
            "",
            "Interpretive cautions:",
            "  - No propagated pointwise measurement uncertainty was supplied.",
            "  - Confidence intervals are conditional on the selected profile and "
            "linear-background model.",
            "  - The residual bootstrap uses short circular blocks when lag-1 "
            "correlation is material; it does not model every possible systematic "
            "error.",
            "  - Apparent FWHM is not instrument-deconvolved and should not be "
            "converted to crystallite size/strain without a resolution standard.",
            "  - 'not_detected' rows retain a fixed-shape area estimate/upper bound; "
            "their q0, d, and FWHM are intentionally not reported.",
            "  - 'detected_fixed_shape' rows have evidence for non-zero area, but "
            "their position/FWHM were borrowed from the configured anchor frame.",
            "  - A 1D radial series cannot establish face-on/edge-on orientation "
            "or mosaicity without the corresponding 2D/azimuthal information.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(args: argparse.Namespace) -> dict[str, Path]:
    input_path = Path(args.input).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.outdir).resolve()
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    config_snapshot_path = output_dir / "config_used.json"
    shutil.copyfile(config_path, config_snapshot_path)

    config = load_config(config_path)
    settings = config["fit_settings"]
    for window in config["windows"]:
        window.setdefault(
            "background_order", int(settings.get("background_order", 1))
        )
    q_all, frames_all, temperatures, _ = load_data(
        input_path,
        {str(key): value for key, value in config["frame_temperatures_C"].items()},
    )
    if args.stage == "first":
        first_frame = (
            canonical_frame(args.first_frame)
            if args.first_frame is not None
            else min(frames_all, key=lambda frame: temperatures[frame])
        )
        if first_frame not in frames_all:
            raise ValueError(
                f"Requested first-frame checkpoint {first_frame} is absent "
                f"from the input."
            )
        frames_all = {first_frame: frames_all[first_frame]}
        temperatures = {first_frame: temperatures[first_frame]}
    elif args.stage == "selected":
        missing_selected = [
            frame for frame in args.frames if frame not in frames_all
        ]
        if missing_selected:
            raise ValueError(
                "Requested representative frames are absent from the input: "
                + ", ".join(missing_selected)
            )
        frames_all = {
            frame: frames_all[frame] for frame in args.frames
        }
        temperatures = {
            frame: temperatures[frame] for frame in args.frames
        }
    dq = float(np.median(np.diff(q_all)))
    plot_series_overview(
        q_all,
        frames_all,
        temperatures,
        figures_dir / "raw_series_overview.png",
    )

    model_selection_rows: list[dict[str, Any]] = []
    parameter_rows: list[dict[str, Any]] = []
    quality_rows: list[dict[str, Any]] = []
    long_rows: list[dict[str, Any]] = []
    final_fits_by_window: dict[str, dict[str, CurveFit]] = {}
    detections_by_window: dict[
        str, dict[str, dict[str, dict[str, Any]]]
    ] = {}

    for window_index, window in enumerate(config["windows"]):
        mask = (q_all >= float(window["q_min"])) & (
            q_all <= float(window["q_max"])
        )
        q = q_all[mask]
        temperature_min = float(window.get("temperature_min_C", -np.inf))
        temperature_max = float(window.get("temperature_max_C", np.inf))
        # Optional explicit acquisition-frame gate (zero-based, inclusive).
        # It complements the temperature gate for isothermal holds in which
        # the pattern changes with time rather than temperature; a window
        # without frame_min/frame_max behaves exactly as before.
        frame_min = int(window.get("frame_min", -1))
        frame_max = int(window.get("frame_max", 10**9))
        window_frames = [
            frame
            for frame in frames_all
            if temperature_min <= temperatures[frame] <= temperature_max
            and frame_min <= int(canonical_frame(frame)) <= frame_max
        ]
        if not window_frames:
            if args.stage == "first":
                print(
                    f"[skip] {window['key']}: checkpoint frame is outside "
                    "the configured temperature/frame range"
                )
                continue
            raise ValueError(
                f"{window['key']}: no frames fall within configured "
                f"temperature range {temperature_min:g} to "
                f"{temperature_max:g} C and frame range "
                f"{frame_min} to {frame_max}"
            )
        frame_data = {
            frame: frames_all[frame][mask] for frame in window_frames
        }
        window_temperatures = {
            frame: temperatures[frame] for frame in window_frames
        }
        if len(q) < 8:
            raise ValueError(
                f"{window['key']}: fitting range contains only {len(q)} q-points"
            )

        (
            selected_profile,
            selected_eta,
            preliminary_fits,
            candidate_rows,
            _,
        ) = compare_profiles(
            q,
            frame_data,
            window,
            window.get(
                "profiles_to_compare", config["profiles_to_compare"]
            ),
            settings,
            seed=args.seed + 10000 * window_index,
            eta_override=(
                float(window["eta_fixed"])
                if "eta_fixed" in window
                else None
            ),
        )
        for row in candidate_rows:
            row["n_frames_fitted"] = len(window_frames)
            row["temperature_min_fitted_C"] = min(
                window_temperatures.values()
            )
            row["temperature_max_fitted_C"] = max(
                window_temperatures.values()
            )
        model_selection_rows.extend(candidate_rows)

        preliminary_detections: dict[
            str, dict[str, dict[str, Any]]
        ] = {}
        for frame, y in frame_data.items():
            preliminary_detections[frame] = test_peak_detection(
                q,
                y,
                window,
                preliminary_fits[frame],
                selected_profile,
                selected_eta,
                settings,
                seed=args.seed
                + 10000 * window_index
                + int(canonical_frame(frame)),
            )
        references = choose_reference_shapes(
            window,
            preliminary_fits,
            preliminary_detections,
            window_temperatures,
        )
        references = refine_reference_shapes_at_anchors(
            q,
            frame_data,
            window,
            preliminary_fits,
            preliminary_detections,
            references,
            window_temperatures,
            selected_profile,
            selected_eta,
            settings,
            seed=args.seed + 50000 * window_index,
        )

        final_fits: dict[str, CurveFit] = {}
        for frame_index, (frame, y) in enumerate(frame_data.items()):
            fixed_shapes = {
                key: {
                    "center": reference["center"],
                    "fwhm": reference["fwhm"],
                }
                for key, reference in references.items()
                if preliminary_detections[frame][key]["status"]
                in {"not_detected", "unresolved"}
                and frame != str(reference["reference_frame"])
            }
            local_seed = (
                args.seed
                + 100000 * window_index
                + 1000 * frame_index
                + int(canonical_frame(frame))
            )
            fit = fit_curve(
                q,
                y,
                window,
                selected_profile,
                eta_fixed=selected_eta,
                fixed_shapes=fixed_shapes,
                initial=preliminary_fits[frame].values,
                n_starts=int(settings["multistarts_final"]),
                rng=np.random.default_rng(local_seed),
                noise_override=preliminary_fits[frame].noise_sigma,
            )
            final_fits[frame] = fit
            for peak in window["peaks"]:
                refine_detection_after_stabilization(
                    preliminary_detections[frame][peak["key"]],
                    fit,
                    peak,
                    frame
                    == str(
                        references[peak["key"]]["reference_frame"]
                    ),
                    settings,
                )
            (
                bootstrap_samples,
                bootstrap_successes,
                bootstrap_block_length,
            ) = bootstrap_fit(
                fit,
                window,
                selected_profile,
                selected_eta,
                fixed_shapes,
                count=args.bootstrap,
                seed=local_seed + 777,
            )
            parameter_rows.extend(
                parameter_rows_for_fit(
                    frame,
                    window_temperatures[frame],
                    window,
                    fit,
                    preliminary_detections[frame],
                    references,
                    bootstrap_samples,
                    bootstrap_successes,
                    args.bootstrap,
                    dq,
                    settings,
                )
            )
            quality_rows.append(
                {
                    **fit_quality_row(
                        frame, window_temperatures[frame], window, fit
                    ),
                    "bootstrap_block_length": bootstrap_block_length,
                }
            )
            long_rows.extend(
                long_rows_for_fit(
                    frame, window_temperatures[frame], window, fit
                )
            )

        final_fits_by_window[window["key"]] = final_fits
        detections_by_window[window["key"]] = preliminary_detections
        plot_window_diagnostics(
            window,
            final_fits,
            preliminary_detections,
            window_temperatures,
            figures_dir / f"fits_{window['key']}.png",
        )
        selected_text = PROFILE_LABELS[selected_profile]
        if selected_profile == "pseudo_voigt":
            selected_text += f" (shared eta={selected_eta:.3f})"
        print(f"[{window_index + 1}/{len(config['windows'])}] {window['key']}: {selected_text}")

    model_selection = pd.DataFrame(model_selection_rows)
    parameters = pd.DataFrame(parameter_rows)
    quality = pd.DataFrame(quality_rows)
    fits_long = pd.DataFrame(long_rows)

    model_selection_path = output_dir / "model_selection.csv"
    parameters_path = output_dir / "peak_parameters.csv"
    quality_path = output_dir / "fit_quality.csv"
    fits_long_path = output_dir / "fitted_curves.csv"
    model_selection.to_csv(model_selection_path, index=False)
    parameters.to_csv(parameters_path, index=False)
    quality.to_csv(quality_path, index=False)
    fits_long.to_csv(fits_long_path, index=False)

    if args.stage == "all":
        plot_parameter_small_multiples(
            parameters,
            config,
            figures_dir / "peak_area_vs_temperature.png",
            metric="area",
        )
        plot_parameter_small_multiples(
            parameters,
            config,
            figures_dir / "peak_position_vs_temperature.png",
            metric="position",
        )
        plot_parameter_small_multiples(
            parameters,
            config,
            figures_dir / "apparent_fwhm_vs_temperature.png",
            metric="fwhm",
        )
        plot_parameter_small_multiples(
            parameters,
            config,
            figures_dir / "peak_area_vs_frame.png",
            metric="area",
            x_axis="frame",
        )
        plot_parameter_small_multiples(
            parameters,
            config,
            figures_dir / "peak_position_vs_frame.png",
            metric="position",
            x_axis="frame",
        )
        plot_parameter_small_multiples(
            parameters,
            config,
            figures_dir / "apparent_fwhm_vs_frame.png",
            metric="fwhm",
            x_axis="frame",
        )
    plot_model_selection(
        model_selection, figures_dir / "profile_model_selection.png"
    )
    plot_qc_heatmap(quality, figures_dir / "fit_quality_heatmap.png")

    checkpoint_path: Path | None = None
    checkpoint_notes_path: Path | None = None
    if args.stage == "first":
        checkpoint_frame = next(iter(frames_all))
        checkpoint_path = figures_dir / "first_frame_diagnostic.png"
        plot_first_frame_diagnostic(
            config,
            final_fits_by_window,
            parameters,
            quality,
            checkpoint_frame,
            temperatures[checkpoint_frame],
            checkpoint_path,
        )
        checkpoint_notes_path = output_dir / "FIRST_FRAME_REVIEW.txt"
        successful = int(quality["optimizer_success"].sum())
        flagged = int(quality["quality_flags"].fillna("").ne("").sum())
        selected_rows = model_selection[model_selection["selected"]]
        checkpoint_lines = [
            "FIRST-FRAME MANUAL REVIEW CHECKPOINT",
            "====================================",
            "",
            f"Frame: {checkpoint_frame}",
            f"Temperature: {temperatures[checkpoint_frame]:g} C",
            f"Optimizer successes: {successful} / {len(quality)} windows",
            f"Windows carrying a QC flag: {flagged} / {len(quality)}",
            "",
            "Selected profile by q-window (Gaussian and Lorentzian were fitted "
            "as alternatives to constrained pseudo-Voigt):",
        ]
        for result in selected_rows.itertuples(index=False):
            eta = (
                f", eta={result.eta_shared:.3f}"
                if result.profile == "pseudo_voigt"
                else ""
            )
            checkpoint_lines.append(
                f"  - {result.window}: {PROFILE_LABELS[result.profile]}{eta}"
            )
        checkpoint_lines.extend(
            [
                "",
                "Manual gate:",
                "  1. Inspect figures/first_frame_diagnostic.png.",
                "  2. Check peak count, bounds, background, component overlap, "
                "and residual structure.",
                "  3. Only after approval, run the all-frame command with "
                "--stage all --confirm-first-fit.",
                "",
                "This stage intentionally does not launch the remaining "
                "frames.",
            ]
        )
        checkpoint_notes_path.write_text(
            "\n".join(checkpoint_lines) + "\n",
            encoding="utf-8",
        )

    summary_path = output_dir / "run_summary.txt"
    write_summary(
        summary_path,
        model_selection,
        parameters,
        quality,
        dq,
        args.bootstrap,
        float(settings["minimum_points_across_fwhm_for_reporting"]),
    )
    metadata = {
        "input_file": str(input_path),
        "config_file": str(config_path),
        "config_snapshot": str(config_snapshot_path),
        "q_step_Ainv": dq,
        "frames": list(frames_all),
        "temperatures_C": temperatures,
        "stage": args.stage,
        "first_frame_confirmation_required_for_all_stage": True,
        "bootstrap_replicates_per_window_frame": args.bootstrap,
        "random_seed": args.seed,
        "python": sys.version,
        "platform": platform.platform(),
        "package_versions": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "method": {
            "input_preprocessing": "none; raw supplied values, no smoothing or offset",
            "background": (
                "simultaneous polynomial in each local q-window; order is "
                "taken from the window configuration"
            ),
            "profile_candidates": config["profiles_to_compare"],
            "selection": (
                "BIC within each window's configured temperature subset; "
                "simpler model within delta BIC < 2"
            ),
            "pseudo_voigt_constraint": "one eta shared across all frames within a q-window",
            "detection": "leave-one-peak-out delta BIC plus area signal-to-noise",
            "uncertainty": (
                "adaptive circular moving-block residual bootstrap when "
                "requested (block length 1 for weak residual correlation; "
                "otherwise approximately n^(1/3)); otherwise local covariance "
                "approximation"
            ),
        },
    }
    metadata_path = output_dir / "run_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    outputs = {
        "config_snapshot": config_snapshot_path,
        "parameters": parameters_path,
        "quality": quality_path,
        "model_selection": model_selection_path,
        "fitted_curves": fits_long_path,
        "summary": summary_path,
        "metadata": metadata_path,
        "figures": figures_dir,
    }
    if checkpoint_path is not None:
        outputs["first_frame_diagnostic"] = checkpoint_path
    if checkpoint_notes_path is not None:
        outputs["first_frame_review"] = checkpoint_notes_path
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Gaussian/Lorentzian/pseudo-Voigt profiles and fit a "
            "multi-temperature GIWAXS line-cut series."
        )
    )
    parser.add_argument("input", help="Tab-delimited text/CSV input (or Excel if supported)")
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
        help="JSON file defining frame temperatures, windows, peaks, and bounds",
    )
    parser.add_argument(
        "--outdir",
        default=str(PROJECT_ROOT / "results" / "drop40" / "example_run"),
        help="Output directory (created if needed)",
    )
    parser.add_argument(
        "--stage",
        choices=("first", "selected", "all"),
        default="first",
        help=(
            "Run the first-frame checkpoint (default), a representative "
            "frame subset, or all configured frames after confirmation."
        ),
    )
    parser.add_argument(
        "--first-frame",
        default=None,
        help=(
            "Frame number for the checkpoint stage. Defaults to the "
            "lowest-temperature frame."
        ),
    )
    parser.add_argument(
        "--frames",
        type=parse_frame_list,
        default=None,
        help=(
            "Comma-separated zero-based frames used with --stage selected, "
            "for example 0,25,50,57,71,85,99."
        ),
    )
    parser.add_argument(
        "--confirm-first-fit",
        action="store_true",
        help=(
            "Required with --stage all; records that the first-frame "
            "diagnostic was manually reviewed."
        ),
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help=(
            "Adaptive moving-block residual-bootstrap replicates per "
            "window/frame. Use 100-500 for final dissertation estimates; "
            "0 uses covariance intervals."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260723,
        help="Random seed for multistart fitting and bootstrap",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.bootstrap < 0:
        parser.error("--bootstrap must be zero or positive")
    if args.stage == "selected" and not args.frames:
        parser.error("--stage selected requires --frames")
    if args.stage != "selected" and args.frames:
        parser.error("--frames is only valid with --stage selected")
    if args.stage == "all" and not args.confirm_first_fit:
        parser.error(
            "--stage all is locked until the diagnostic is reviewed; rerun "
            "with --confirm-first-fit after approval."
        )
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(args.outdir).resolve() / ".matplotlib")
    )
    outputs = run_pipeline(args)
    print("\nCompleted peak fitting.")
    for label, path in outputs.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
