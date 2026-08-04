#!/usr/bin/env python3
"""Build scan-specific candidate-v2 fitting configs for GIWAXS scan 587225.

The windows below come from visual review of the unsmoothed, d5i-normalized
FR/IP/OOP profiles at frames 0, 10, 19, 21, 24, 27, 30, 40, 70, and 134.
They are intentionally not copied from scan 587214.  Temperature gates separate
features observed through frame 27 (63.9 C) from those established by frame 30
(71.8 C), without asserting a thermodynamic transition temperature.

This script writes configuration files only.  Candidate runs must use
``--bootstrap 0`` and pass visual first-frame checkpoints before selected-frame
fitting.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = (
    PROJECT_ROOT
    / "data"
    / "dimitar_n2n2"
    / "scan_587225"
    / "processed"
    / "scan587225_frame_manifest.csv"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT
    / "configs"
    / "secondary_samples"
    / "scan_587225"
    / "candidate_v2"
)
SELECTED_FRAMES = [0, 10, 19, 21, 24, 27, 30, 40, 70, 134]
PROFILE_FAMILIES = ["gaussian", "lorentzian", "pseudo_voigt"]
EARLY_MAX_C = 64.9
LATE_MIN_C = 65.0


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
    anchor_temperature_C: float,
) -> dict[str, Any]:
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
        "anchor_temperature_C": anchor_temperature_C,
        "provisional_assignment": True,
    }


def window(
    key: str,
    title: str,
    q_min: float,
    q_max: float,
    peaks: list[dict[str, Any]],
    rationale: str,
    *,
    temperature_min_C: float | None = None,
    temperature_max_C: float | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "key": key,
        "title": title,
        "q_min": q_min,
        "q_max": q_max,
        "peaks": peaks,
        "profiles_to_compare": PROFILE_FAMILIES,
        "candidate_window_requires_visual_review": True,
        "stage1_rationale": rationale,
    }
    if temperature_min_C is not None:
        result["temperature_min_C"] = temperature_min_C
    if temperature_max_C is not None:
        result["temperature_max_C"] = temperature_max_C
    return result


def common_low_q(cut: str) -> dict[str, Any]:
    return window(
        "low_q_pair",
        f"{cut} provisional low-q broad/0.331 pair",
        0.225,
        0.405,
        [
            peak(
                "low_q_broad",
                f"{cut} broad low-q component near 0.30",
                0.305,
                0.285,
                0.321,
                0.075,
                0.035,
                0.130,
                40,
                84.1,
            ),
            peak(
                "low_q_0331",
                f"{cut} sharper low-q component near 0.331",
                0.331,
                0.324,
                0.340,
                0.018,
                0.009,
                0.045,
                40,
                84.1,
            ),
        ],
        (
            "All selected profiles show a broad envelope near 0.30 and a local "
            "maximum near 0.331 A^-1.  The inherited 0.338-0.360 A^-1 bound "
            "would miss this scan's sharper maximum."
        ),
    )


def fr_windows() -> list[dict[str, Any]]:
    return [
        common_low_q("FR"),
        window(
            "mid_q_early_pair",
            "FR low-temperature/ramp features near 0.436 and 0.486",
            0.405,
            0.515,
            [
                peak(
                    "feature_0436_early",
                    "FR early feature near 0.436",
                    0.436,
                    0.423,
                    0.450,
                    0.020,
                    0.009,
                    0.050,
                    0,
                    31.3,
                ),
                peak(
                    "feature_0486_early_ramp",
                    "FR weak early/ramp shoulder near 0.486",
                    0.486,
                    0.468,
                    0.501,
                    0.035,
                    0.015,
                    0.080,
                    0,
                    31.3,
                ),
            ],
            (
                "The 0.436 A^-1 peak is clear at low temperature; a repeatable "
                "weak shoulder occurs near 0.486 A^-1.  The latter remains "
                "optional and must pass residual/evidence review."
            ),
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "mid_q_late_pair",
            "FR late features near 0.455 and 0.533",
            0.425,
            0.580,
            [
                peak(
                    "feature_0455_late",
                    "FR late feature near 0.455",
                    0.455,
                    0.443,
                    0.467,
                    0.022,
                    0.009,
                    0.055,
                    40,
                    84.1,
                ),
                peak(
                    "feature_0533_late",
                    "FR late broad feature near 0.533",
                    0.533,
                    0.515,
                    0.550,
                    0.045,
                    0.018,
                    0.100,
                    40,
                    84.1,
                ),
            ],
            (
                "By frame 30/71.8 C the early pair is replaced by resolved "
                "structure near 0.455 and 0.530-0.536 A^-1."
            ),
            temperature_min_C=LATE_MIN_C,
        ),
        window(
            "feature_0754_early",
            "FR early peak near 0.754",
            0.700,
            0.810,
            [
                peak(
                    "feature_0754_early",
                    "FR early peak near 0.754",
                    0.754,
                    0.742,
                    0.765,
                    0.018,
                    0.009,
                    0.050,
                    0,
                    31.3,
                )
            ],
            "A sharp 0.754 A^-1 peak is visible through the early ramp and then disappears.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_0941_early_broad",
            "FR weak early broad feature near 0.941",
            0.870,
            1.020,
            [
                peak(
                    "feature_0941_early_broad",
                    "FR weak broad feature near 0.941",
                    0.941,
                    0.915,
                    0.965,
                    0.065,
                    0.025,
                    0.120,
                    0,
                    31.3,
                )
            ],
            "The low-temperature FR cut has a broad, weak maximum near 0.941 A^-1.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_1535_early",
            "FR early peak near 1.535",
            1.470,
            1.590,
            [
                peak(
                    "feature_1535_early",
                    "FR early peak near 1.532-1.538",
                    1.535,
                    1.518,
                    1.550,
                    0.020,
                    0.009,
                    0.055,
                    0,
                    31.3,
                )
            ],
            "A sharp 1.532-1.538 A^-1 peak weakens to the noise floor before the high-temperature period.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_1720_early",
            "FR early broad feature near 1.72",
            1.650,
            1.800,
            [
                peak(
                    "feature_1720_early",
                    "FR early feature near 1.71-1.73",
                    1.720,
                    1.695,
                    1.765,
                    0.060,
                    0.018,
                    0.120,
                    0,
                    31.3,
                )
            ],
            "A broad/evolving early feature lies above the inherited 1.700 A^-1 centre bounds.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "high_q_early_pair",
            "FR early peaks near 1.999 and 2.08",
            1.930,
            2.160,
            [
                peak(
                    "feature_1999_early",
                    "FR early peak near 1.999",
                    1.999,
                    1.980,
                    2.015,
                    0.022,
                    0.009,
                    0.055,
                    0,
                    31.3,
                ),
                peak(
                    "feature_2080_early",
                    "FR early feature near 2.08",
                    2.080,
                    2.050,
                    2.105,
                    0.045,
                    0.015,
                    0.090,
                    0,
                    31.3,
                ),
            ],
            "Two separated high-q features are visible at low temperature; the inherited single 2.005 A^-1 component is insufficient.",
            temperature_max_C=EARLY_MAX_C,
        ),
    ]


def ip_windows() -> list[dict[str, Any]]:
    return [
        common_low_q("IP"),
        window(
            "mid_q_moving",
            "IP evolving feature from about 0.486 to 0.455",
            0.420,
            0.530,
            [
                peak(
                    "feature_0486_to_0455",
                    "IP feature moving from about 0.486 to 0.455",
                    0.455,
                    0.445,
                    0.495,
                    0.030,
                    0.012,
                    0.080,
                    40,
                    84.1,
                )
            ],
            "One dominant maximum moves from about 0.486 A^-1 at low temperature to about 0.455 A^-1 at high temperature; start with one component rather than forcing the inherited pair.",
        ),
        window(
            "feature_1544_early",
            "IP early peak near 1.544",
            1.470,
            1.590,
            [
                peak(
                    "feature_1544_early",
                    "IP early peak near 1.544",
                    1.544,
                    1.525,
                    1.560,
                    0.020,
                    0.009,
                    0.060,
                    0,
                    31.3,
                )
            ],
            "Only one distinct low-temperature peak is visible in this range; no separate 1.495 A^-1 component is justified.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_1720_early",
            "IP early peak near 1.72",
            1.660,
            1.790,
            [
                peak(
                    "feature_1720_early",
                    "IP early feature near 1.71-1.73",
                    1.720,
                    1.695,
                    1.742,
                    0.040,
                    0.012,
                    0.085,
                    0,
                    31.3,
                )
            ],
            "The early maximum occurs below the inherited 1.744 A^-1 centre bounds.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_2080_early",
            "IP early broad peak near 2.08",
            1.980,
            2.150,
            [
                peak(
                    "feature_2080_early",
                    "IP early broad peak near 2.08",
                    2.080,
                    2.045,
                    2.105,
                    0.060,
                    0.018,
                    0.120,
                    0,
                    31.3,
                )
            ],
            "The dominant IP high-q peak is near 2.08 A^-1 and lies outside the inherited 2.011 A^-1 window.",
            temperature_max_C=EARLY_MAX_C,
        ),
    ]


def oop_windows() -> list[dict[str, Any]]:
    return [
        common_low_q("OOP"),
        window(
            "feature_0436_early",
            "OOP early peak near 0.436",
            0.405,
            0.485,
            [
                peak(
                    "feature_0436_early",
                    "OOP early peak near 0.436",
                    0.436,
                    0.423,
                    0.450,
                    0.020,
                    0.009,
                    0.050,
                    0,
                    31.3,
                )
            ],
            "A strong low-temperature OOP peak near 0.436 A^-1 fades during the ramp.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_0536_late",
            "OOP late broad peak near 0.536",
            0.490,
            0.585,
            [
                peak(
                    "feature_0536_late",
                    "OOP late broad peak near 0.536",
                    0.536,
                    0.515,
                    0.552,
                    0.045,
                    0.018,
                    0.100,
                    40,
                    84.1,
                )
            ],
            "A separate broad OOP feature near 0.536 A^-1 is established after the ramp.",
            temperature_min_C=LATE_MIN_C,
        ),
        window(
            "feature_0754_early",
            "OOP early peak near 0.754",
            0.700,
            0.800,
            [
                peak(
                    "feature_0754_early",
                    "OOP early peak near 0.754",
                    0.754,
                    0.742,
                    0.766,
                    0.018,
                    0.009,
                    0.050,
                    0,
                    31.3,
                )
            ],
            "The strongest early OOP peak is at about 0.754 A^-1 and disappears by the high-temperature period.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_0835",
            "OOP feature near 0.835",
            0.790,
            0.875,
            [
                peak(
                    "feature_0835",
                    "OOP feature near 0.829-0.841",
                    0.835,
                    0.815,
                    0.850,
                    0.030,
                    0.012,
                    0.070,
                    40,
                    84.1,
                )
            ],
            "A weaker feature near 0.835 A^-1 is visible both early and at frame 40; it is distinct from the 0.754 A^-1 peak.",
        ),
        window(
            "feature_1495_early",
            "OOP early peak near 1.495",
            1.440,
            1.570,
            [
                peak(
                    "feature_1495_early",
                    "OOP early peak near 1.495",
                    1.495,
                    1.478,
                    1.512,
                    0.020,
                    0.009,
                    0.060,
                    0,
                    31.3,
                )
            ],
            "A sharp low-temperature OOP peak near 1.495 A^-1 is absent from the generic OOP template.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "high_q_early_pair",
            "OOP early peaks near 1.71 and 1.75",
            1.650,
            1.800,
            [
                peak(
                    "feature_1710_early",
                    "OOP early peak near 1.71",
                    1.710,
                    1.690,
                    1.726,
                    0.022,
                    0.009,
                    0.060,
                    0,
                    31.3,
                ),
                peak(
                    "feature_1750_early",
                    "OOP secondary early peak near 1.75",
                    1.750,
                    1.735,
                    1.772,
                    0.025,
                    0.009,
                    0.075,
                    0,
                    31.3,
                ),
            ],
            "The first OOP frame resolves maxima near 1.71 and 1.75 A^-1; disjoint centre bounds prevent swapping.",
            temperature_max_C=EARLY_MAX_C,
        ),
        window(
            "feature_2223_early",
            "OOP early peak near 2.223",
            2.150,
            2.290,
            [
                peak(
                    "feature_2223_early",
                    "OOP early peak near 2.223",
                    2.223,
                    2.198,
                    2.240,
                    0.026,
                    0.009,
                    0.070,
                    0,
                    31.3,
                )
            ],
            "A distinct low-temperature OOP peak near 2.223 A^-1 is absent from the generic template.",
            temperature_max_C=EARLY_MAX_C,
        ),
    ]


def validate_windows(windows: list[dict[str, Any]]) -> None:
    seen_windows: set[str] = set()
    seen_peaks: set[str] = set()
    for item in windows:
        if item["key"] in seen_windows:
            raise ValueError(f"Duplicate window key: {item['key']}")
        seen_windows.add(item["key"])
        if item["q_min"] >= item["q_max"]:
            raise ValueError(f"Invalid q window: {item['key']}")
        intervals = []
        for component in item["peaks"]:
            key = component["key"]
            if key in seen_peaks:
                raise ValueError(f"Duplicate peak key within one cut: {key}")
            seen_peaks.add(key)
            if not (
                item["q_min"]
                < component["q_min"]
                < component["q_guess"]
                < component["q_max"]
                < item["q_max"]
            ):
                raise ValueError(f"Peak bounds escape window for {key}")
            if not (
                0
                < component["fwhm_min"]
                <= component["fwhm_guess"]
                <= component["fwhm_max"]
            ):
                raise ValueError(f"Invalid FWHM bounds for {key}")
            intervals.append((component["q_min"], component["q_max"], key))
        intervals.sort()
        for left, right in zip(intervals, intervals[1:]):
            if left[1] >= right[0]:
                raise ValueError(
                    f"Overlapping centre bounds in {item['key']}: "
                    f"{left[2]} and {right[2]}"
                )


def build_config(
    cut: str,
    windows: list[dict[str, Any]],
    manifest: pd.DataFrame,
    manifest_path: Path,
) -> dict[str, Any]:
    validate_windows(windows)
    temperatures = {
        str(int(row.frame)): float(row.measured_temperature_C)
        for row in manifest.itertuples(index=False)
    }
    return {
        "frame_temperatures_C": temperatures,
        "profiles_to_compare": PROFILE_FAMILIES,
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
        "analysis_scope": (
            "Scan-specific candidate-v2 key-frame screening for peak presence, "
            "q0/d and integrated area. No CCL interpretation."
        ),
        "representative_frames": SELECTED_FRAMES,
        "windows": windows,
        "series_metadata": {
            "scan": 587225,
            "sample_label": "Dimitar_n2n2",
            "cut": cut,
            "candidate_version": 2,
            "series_type": "secondary_in_situ_time_series",
            "frame_manifest": str(manifest_path),
            "key_frames": SELECTED_FRAMES,
            "stage1_plot_root": str(
                PROJECT_ROOT
                / "results"
                / "dimitar_n2n2"
                / "scan_587225"
                / "key_frames_v1"
            ),
            "source_of_windows": (
                "Independent visual review of unsmoothed d5i-normalized scan "
                "587225 profiles; not a transfer of scan 587214 windows."
            ),
            "temperature_gate_note": (
                "64.9/65.0 C gates separate the last selected frame with the "
                "early pattern (frame 27, 63.9 C) from the first selected frame "
                "with the late pattern (frame 30, 71.8 C). They are computational "
                "gates, not claimed phase-transition temperatures."
            ),
            "bootstrap_policy": "Candidate and checkpoint runs use --bootstrap 0.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build visually reviewed candidate-v2 configs for scan 587225."
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    manifest_path = args.manifest.expanduser().resolve()
    manifest = pd.read_csv(manifest_path)
    required = {"frame", "measured_temperature_C", "key_frame_selected"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise ValueError(f"Manifest is missing: {', '.join(missing)}")
    manifest["frame"] = pd.to_numeric(manifest["frame"], errors="raise").astype(int)
    manifest["measured_temperature_C"] = pd.to_numeric(
        manifest["measured_temperature_C"], errors="raise"
    ).astype(float)
    if manifest["frame"].tolist() != list(range(135)):
        raise ValueError("Expected continuous frames 0-134 for scan 587225.")
    selected = manifest.loc[
        manifest["key_frame_selected"].astype(str).str.lower().isin(
            {"true", "1", "yes"}
        ),
        "frame",
    ].tolist()
    if selected != SELECTED_FRAMES:
        raise ValueError(
            f"Manifest key frames {selected} do not match {SELECTED_FRAMES}."
        )

    configs = {
        "FR": build_config("FR", fr_windows(), manifest, manifest_path),
        "IP": build_config("IP", ip_windows(), manifest, manifest_path),
        "OOP": build_config("OOP", oop_windows(), manifest, manifest_path),
    }
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for cut, config in configs.items():
        output = output_dir / f"scan587225_{cut}_candidate_keyframes.json"
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
