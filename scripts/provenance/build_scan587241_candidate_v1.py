#!/usr/bin/env python3
"""Build the auditable scan-specific candidate-v1 configs for scan 587241.

The candidate windows are transcribed from the saved Stage 1 visual review.
This script writes configuration only; it does not run peak fitting.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCAN = "587241"
FRAMES = [0, 9, 18, 20, 25, 30, 35, 40, 45, 50, 60, 80, 129]
EARLY_ANCHOR = 18
TRANSFORMED_ANCHOR = 40
MANIFEST = (
    PROJECT_ROOT
    / "data"
    / "dimitar_airn2"
    / "scan_587241"
    / "processed"
    / "scan587241_frame_manifest.csv"
)
VISUAL_REVIEW = (
    PROJECT_ROOT
    / "results"
    / "dimitar_airn2"
    / "scan_587241"
    / "key_frames_v1"
    / "STAGE1_VISUAL_REVIEW.md"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "configs"
    / "secondary_samples"
    / "scan_587241"
    / "candidate_v1"
)


def peak(
    key: str,
    label: str,
    q_guess: float,
    q_min: float,
    q_max: float,
    fwhm_guess: float,
    fwhm_min: float,
    fwhm_max: float,
    anchor_frame: int,
    temperatures: dict[str, float],
) -> dict[str, object]:
    return {
        "key": key,
        "label": label,
        "q_guess": q_guess,
        "q_min": q_min,
        "q_max": q_max,
        "fwhm_guess": fwhm_guess,
        "fwhm_min": fwhm_min,
        "fwhm_max": fwhm_max,
        "anchor_frame": anchor_frame,
        "anchor_temperature_C": temperatures[str(anchor_frame)],
        "provisional_assignment": True,
    }


def window(
    key: str,
    title: str,
    q_min: float,
    q_max: float,
    peaks: list[dict[str, object]],
    regime: str,
    *,
    temperature_min_C: float | None = None,
    temperature_max_C: float | None = None,
) -> dict[str, object]:
    result: dict[str, object] = {
        "key": key,
        "title": title,
        "q_min": q_min,
        "q_max": q_max,
        "background_order": 1,
        "peaks": peaks,
        "regime_label": regime,
        "candidate_window_requires_visual_review": True,
    }
    if temperature_min_C is not None:
        result["temperature_min_C"] = temperature_min_C
    if temperature_max_C is not None:
        result["temperature_max_C"] = temperature_max_C
    return result


def fr_windows(temperatures: dict[str, float]) -> list[dict[str, object]]:
    early = EARLY_ANCHOR
    late = TRANSFORMED_ANCHOR
    return [
        window(
            "fr_early_q0331",
            "FR early low-q feature near 0.331",
            0.275,
            0.375,
            [peak("fr_q0331", "FR unassigned early component near 0.331", 0.331, 0.312, 0.344, 0.025, 0.009, 0.070, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "fr_early_q0430",
            "FR early sharp feature near 0.430",
            0.390,
            0.470,
            [peak("fr_q0430", "FR unassigned early sharp component near 0.430", 0.430, 0.416, 0.444, 0.015, 0.009, 0.040, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "fr_early_q0742",
            "FR early sharp feature near 0.742",
            0.690,
            0.795,
            [peak("fr_q0742", "FR unassigned early component near 0.742", 0.742, 0.725, 0.760, 0.018, 0.009, 0.050, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "fr_early_q1510",
            "FR early high-q feature near 1.51",
            1.445,
            1.575,
            [peak("fr_q1510", "FR unassigned early component near 1.51", 1.510, 1.485, 1.545, 0.032, 0.012, 0.080, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "fr_early_q1695",
            "FR early high-q feature near 1.70",
            1.635,
            1.760,
            [peak("fr_q1695", "FR unassigned early component near 1.70", 1.695, 1.670, 1.720, 0.030, 0.012, 0.075, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "fr_early_q1965",
            "FR early high-q feature near 1.96",
            1.905,
            2.025,
            [peak("fr_q1965", "FR unassigned early component near 1.96", 1.965, 1.940, 2.005, 0.032, 0.012, 0.080, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "fr_transformed_q0306",
            "FR transformed broad 0.306 feature and 0.331 shoulder",
            0.235,
            0.400,
            [
                peak("fr_q0306", "FR unassigned transformed broad component near 0.306", 0.306, 0.285, 0.322, 0.068, 0.035, 0.115, late, temperatures),
                peak("fr_q0331_shoulder", "FR unassigned transformed shoulder near 0.331", 0.331, 0.324, 0.342, 0.018, 0.009, 0.040, late, temperatures),
            ],
            "transformed_state_ge_60C",
            temperature_min_C=60.0,
        ),
        window(
            "fr_transformed_q0520",
            "FR transformed feature near 0.520",
            0.455,
            0.590,
            [peak("fr_q0520", "FR unassigned transformed component near 0.520", 0.520, 0.495, 0.545, 0.040, 0.015, 0.095, late, temperatures)],
            "transformed_state_ge_60C",
            temperature_min_C=60.0,
        ),
    ]


def ip_windows(temperatures: dict[str, float]) -> list[dict[str, object]]:
    early = EARLY_ANCHOR
    late = TRANSFORMED_ANCHOR
    return [
        window(
            "ip_q0324",
            "IP recurrent low-q feature near 0.324",
            0.275,
            0.380,
            [peak("ip_q0324", "IP unassigned recurrent component near 0.324", 0.324, 0.305, 0.338, 0.020, 0.009, 0.065, late, temperatures)],
            "all_selected_frames",
        ),
        window(
            "ip_early_low_q_pair",
            "IP early low-q components near 0.480 and 0.536",
            0.425,
            0.570,
            [
                peak("ip_q0480", "IP unassigned early low-q component near 0.480", 0.480, 0.463, 0.492, 0.018, 0.009, 0.050, early, temperatures),
                peak("ip_q0536", "IP unassigned early low-q component near 0.536", 0.536, 0.520, 0.550, 0.018, 0.009, 0.050, early, temperatures),
            ],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "ip_transformed_low_q_pair",
            "IP transformed low-q components near 0.449 and 0.505",
            0.405,
            0.565,
            [
                peak("ip_q0449", "IP unassigned transformed low-q component near 0.449", 0.449, 0.438, 0.463, 0.018, 0.009, 0.050, late, temperatures),
                peak("ip_q0505", "IP unassigned transformed low-q component near 0.505", 0.505, 0.491, 0.520, 0.024, 0.009, 0.065, late, temperatures),
            ],
            "transformed_state_ge_60C",
            temperature_min_C=60.0,
        ),
        window(
            "ip_early_q1515",
            "IP early high-q feature near 1.515",
            1.440,
            1.580,
            [peak("ip_q1515", "IP unassigned early component near 1.515", 1.515, 1.485, 1.548, 0.038, 0.012, 0.085, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "ip_early_high_q_pair",
            "IP early high-q shoulder and main feature",
            1.620,
            1.760,
            [
                peak("ip_q1669", "IP unassigned early shoulder near 1.669", 1.669, 1.645, 1.682, 0.035, 0.012, 0.080, early, temperatures),
                peak("ip_q1698", "IP unassigned early main component near 1.698", 1.698, 1.684, 1.720, 0.026, 0.009, 0.060, early, temperatures),
            ],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "ip_early_q1995",
            "IP early high-q feature near 2.00",
            1.940,
            2.025,
            [peak("ip_q1995", "IP unassigned early component near 2.00", 1.995, 1.970, 2.018, 0.028, 0.012, 0.060, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "ip_transformed_q1738",
            "IP transformed weak high-q feature near 1.738",
            1.680,
            1.800,
            [peak("ip_q1738", "IP unassigned transformed component near 1.738", 1.738, 1.715, 1.755, 0.038, 0.015, 0.090, late, temperatures)],
            "transformed_state_ge_60C_weak_candidate",
            temperature_min_C=60.0,
        ),
    ]


def oop_windows(temperatures: dict[str, float]) -> list[dict[str, object]]:
    early = EARLY_ANCHOR
    late = TRANSFORMED_ANCHOR
    return [
        window(
            "oop_q0306",
            "OOP broad low-q feature near 0.306",
            0.235,
            0.390,
            [peak("oop_q0306", "OOP unassigned broad component near 0.306", 0.306, 0.285, 0.325, 0.068, 0.035, 0.115, late, temperatures)],
            "all_selected_frames",
        ),
        window(
            "oop_early_q0430",
            "OOP early sharp feature near 0.430",
            0.390,
            0.470,
            [peak("oop_q0430", "OOP unassigned early sharp component near 0.430", 0.430, 0.416, 0.444, 0.015, 0.009, 0.040, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_early_low_q_pair",
            "OOP early low-q components near 0.55 and 0.59",
            0.490,
            0.625,
            [
                peak("oop_q0548", "OOP unassigned early low-q component A", 0.548, 0.505, 0.565, 0.024, 0.009, 0.080, early, temperatures),
                peak("oop_q0588", "OOP unassigned early low-q component B", 0.588, 0.570, 0.605, 0.020, 0.009, 0.050, early, temperatures),
            ],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_early_q0742",
            "OOP early sharp feature near 0.742",
            0.690,
            0.795,
            [peak("oop_q0742", "OOP unassigned early component near 0.742", 0.742, 0.725, 0.760, 0.018, 0.009, 0.050, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_early_q1478",
            "OOP early high-q feature near 1.478",
            1.410,
            1.540,
            [peak("oop_q1478", "OOP unassigned early component near 1.478", 1.478, 1.445, 1.505, 0.032, 0.012, 0.080, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_early_q1695",
            "OOP early high-q feature near 1.70",
            1.630,
            1.760,
            [peak("oop_q1695", "OOP unassigned early component near 1.70", 1.695, 1.665, 1.720, 0.022, 0.009, 0.070, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_early_q1955",
            "OOP early high-q feature near 1.955",
            1.900,
            2.020,
            [peak("oop_q1955", "OOP unassigned early component near 1.955", 1.955, 1.930, 1.980, 0.032, 0.012, 0.080, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_early_q2198",
            "OOP early high-q feature near 2.198",
            2.140,
            2.250,
            [peak("oop_q2198", "OOP unassigned early component near 2.198", 2.198, 2.170, 2.220, 0.038, 0.012, 0.080, early, temperatures)],
            "early_state_le_55C",
            temperature_max_C=55.0,
        ),
        window(
            "oop_transformed_q0528",
            "OOP transformed feature near 0.528",
            0.460,
            0.600,
            [peak("oop_q0528", "OOP unassigned transformed component near 0.528", 0.528, 0.505, 0.550, 0.038, 0.018, 0.100, late, temperatures)],
            "transformed_state_ge_60C",
            temperature_min_C=60.0,
        ),
    ]


def build_config(
    cut: str,
    temperatures: dict[str, float],
    windows: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "frame_temperatures_C": temperatures,
        "profiles_to_compare": ["gaussian", "lorentzian", "pseudo_voigt"],
        "fit_settings": {
            "background_order": 1,
            "multistarts_model_selection": 4,
            "multistarts_final": 6,
            "pseudo_voigt_eta_grid_step": 0.1,
            "minimum_points_across_fwhm_for_reporting": 5.0,
            "tentative_delta_bic": 2.0,
            "detected_delta_bic": 6.0,
            "strong_delta_bic": 10.0,
            "tentative_area_snr": 2.0,
            "detected_area_snr": 3.0,
            "strong_area_snr": 5.0,
        },
        "windows": windows,
        "analysis_scope": (
            "Sample-specific key-frame screening for peak presence, conservative "
            "q0/d and integrated area. No CCL interpretation; apparent FWHM is "
            "a fitting/QC quantity only."
        ),
        "representative_frames": FRAMES,
        "series_metadata": {
            "scan": int(SCAN),
            "sample_label": "Dimitar_airn2",
            "cut": cut,
            "series_type": "secondary_in_situ_time_series",
            "frame_manifest": str(MANIFEST.resolve()),
            "key_frames": FRAMES,
            "early_regime_temperature_max_C": 55.0,
            "transformed_regime_temperature_min_C": 60.0,
            "early_anchor_frame": EARLY_ANCHOR,
            "early_anchor_temperature_C": temperatures[str(EARLY_ANCHOR)],
            "transformed_anchor_frame": TRANSFORMED_ANCHOR,
            "transformed_anchor_temperature_C": temperatures[str(TRANSFORMED_ANCHOR)],
            "visual_review": str(VISUAL_REVIEW.resolve()),
            "config_builder": str(Path(__file__).resolve()),
            "scientific_status": (
                "Candidate v1 transcribed from the scan-specific Stage 1 raw-profile "
                "review. Early and transformed regimes are separated where peak "
                "identity across the transition is ambiguous. Frame-18/frame-40 "
                "checkpoint residuals motivated the transformed FR 0.331 shoulder, "
                "regime-specific IP low-q pairs, and the early IP 1.669/1.698 pair; "
                "weak optional maxima remain deliberately excluded."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if not MANIFEST.is_file() or not VISUAL_REVIEW.is_file():
        raise FileNotFoundError("Prepared manifest or Stage 1 visual review is missing.")
    manifest = pd.read_csv(MANIFEST)
    if manifest["frame"].astype(int).tolist() != list(range(130)):
        raise ValueError("Scan 587241 manifest is not a continuous 130-frame series.")
    selected = manifest.loc[manifest["key_frame_selected"].astype(bool), "frame"].astype(int).tolist()
    if selected != FRAMES:
        raise ValueError(f"Confirmed frames changed: expected {FRAMES}, found {selected}")
    temperatures = {
        str(int(row.frame)): float(row.measured_temperature_C)
        for row in manifest.itertuples(index=False)
    }
    if temperatures[str(EARLY_ANCHOR)] > 55.0:
        raise ValueError("Early anchor is outside the <=55 C regime.")
    if temperatures[str(TRANSFORMED_ANCHOR)] < 60.0:
        raise ValueError("Transformed anchor is outside the >=60 C regime.")

    builders = {"FR": fr_windows, "IP": ip_windows, "OOP": oop_windows}
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for cut, builder in builders.items():
        output = output_dir / f"scan{SCAN}_{cut}_candidate_keyframes.json"
        if output.exists() and not args.overwrite:
            raise FileExistsError(f"Refusing to overwrite {output}; pass --overwrite after review.")
        config = build_config(cut, temperatures, builder(temperatures))
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
