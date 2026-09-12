#!/usr/bin/env python3
"""Build auditable scan-specific candidate-v2 configs for scan 587217.

The window library is derived from the unsmoothed key-frame visual review in
``results/dimitar_bcnogas1/scan_587217/key_frames_v1``.  It deliberately omits
unsupported scan-587214 template peaks and assigns frame 0 to features visible
only before the rapid heating transition.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCAN = "587217"
SAMPLE_LABEL = "Dimitar_bcnogas1"
KEY_FRAMES = [0, 1, 2, 3, 4, 5, 6, 11, 25, 35]
MANIFEST_PATH = (
    PROJECT_ROOT
    / "data"
    / "dimitar_bcnogas1"
    / "scan_587217"
    / "processed"
    / "scan587217_frame_manifest.csv"
)
VISUAL_REVIEW_PATH = (
    PROJECT_ROOT
    / "results"
    / "dimitar_bcnogas1"
    / "scan_587217"
    / "key_frames_v1"
    / "STAGE1_VISUAL_REVIEW.md"
)
OUTPUT_DIR = (
    PROJECT_ROOT
    / "configs"
    / "secondary_samples"
    / "scan_587217"
    / "candidate_v2"
)

PROFILES = ["gaussian", "lorentzian", "pseudo_voigt"]
FIT_SETTINGS = {
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
}


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
        "anchor_temperature_C": temperatures[str(anchor_frame)],
        "provisional_assignment": True,
    }


def window(
    key: str,
    title: str,
    q_min: float,
    q_max: float,
    peaks: list[dict[str, Any]],
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
        "candidate_window_requires_visual_review": True,
    }
    if temperature_min_C is not None:
        result["temperature_min_C"] = temperature_min_C
    if temperature_max_C is not None:
        result["temperature_max_C"] = temperature_max_C
    return result


def fr_windows(temperatures: dict[str, float]) -> list[dict[str, Any]]:
    return [
        window(
            "low_q_pair",
            "FR broad q=0.306 and sharper q=0.331 components",
            0.245,
            0.375,
            [
                peak("low_q_broad", "FR broad low-q component", 0.306, 0.285, 0.318, 0.065, 0.030, 0.120, 11, temperatures),
                peak("low_q_sharp", "FR sharper low-q component", 0.331, 0.322, 0.340, 0.016, 0.006, 0.035, 11, temperatures),
            ],
        ),
        window(
            "feature_0430_to_0449",
            "FR evolving q=0.430 to 0.449 feature",
            0.395,
            0.480,
            [peak("feature_0430_to_0449", "FR evolving mid-q component", 0.449, 0.416, 0.477, 0.020, 0.006, 0.045, 11, temperatures)],
        ),
        window(
            "feature_0530",
            "FR high-temperature feature near q=0.530",
            0.485,
            0.580,
            [peak("feature_0530", "FR component near q=0.530", 0.530, 0.510, 0.545, 0.040, 0.015, 0.090, 11, temperatures)],
            temperature_min_C=65.0,
        ),
        window(
            "feature_0741_early",
            "FR early feature near q=0.741",
            0.700,
            0.790,
            [peak("feature_0741_early", "FR early component near q=0.741", 0.741, 0.730, 0.753, 0.020, 0.006, 0.040, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_1490_1526_early",
            "FR early q=1.49 shoulder and q=1.526 peak",
            1.470,
            1.580,
            [
                peak("feature_1490_early", "FR early shoulder near q=1.49", 1.490, 1.478, 1.508, 0.035, 0.012, 0.075, 0, temperatures),
                peak("feature_1526_early", "FR early component near q=1.526", 1.526, 1.512, 1.540, 0.024, 0.006, 0.050, 0, temperatures),
            ],
            temperature_max_C=60.0,
        ),
        window(
            "feature_1715_early",
            "FR early feature near q=1.715",
            1.660,
            1.780,
            [peak("feature_1715_early", "FR early component near q=1.715", 1.715, 1.697, 1.730, 0.038, 0.008, 0.070, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_1986_early",
            "FR early feature near q=1.986",
            1.930,
            2.020,
            [peak("feature_1986_early", "FR early sharp component near q=1.986", 1.986, 1.970, 2.002, 0.024, 0.006, 0.050, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_2060_early",
            "FR early broad feature near q=2.060",
            2.025,
            2.120,
            [peak("feature_2060_early", "FR early broad component near q=2.060", 2.060, 2.035, 2.095, 0.055, 0.018, 0.120, 0, temperatures)],
            temperature_max_C=60.0,
        ),
    ]


def ip_windows(temperatures: dict[str, float]) -> list[dict[str, Any]]:
    return [
        window(
            "feature_0327",
            "IP persistent feature near q=0.327",
            0.295,
            0.380,
            [peak("feature_0327", "IP component near q=0.327", 0.327, 0.318, 0.342, 0.018, 0.006, 0.040, 11, temperatures)],
        ),
        window(
            "feature_0474_to_0449",
            "IP evolving q=0.474 to 0.449 feature",
            0.405,
            0.520,
            [peak("feature_0474_to_0449", "IP evolving mid-q component", 0.449, 0.438, 0.482, 0.022, 0.006, 0.055, 11, temperatures)],
        ),
        window(
            "feature_1713_early",
            "IP early feature near q=1.713",
            1.675,
            1.765,
            [peak("feature_1713_early", "IP early component near q=1.713", 1.713, 1.695, 1.725, 0.018, 0.006, 0.050, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_2055_early",
            "IP early broad feature near q=2.055",
            1.980,
            2.130,
            [peak("feature_2055_early", "IP early broad component near q=2.055", 2.055, 2.030, 2.085, 0.065, 0.025, 0.130, 0, temperatures)],
            temperature_max_C=60.0,
        ),
    ]


def oop_windows(temperatures: dict[str, float]) -> list[dict[str, Any]]:
    return [
        window(
            "low_q_broad",
            "OOP persistent broad low-q envelope",
            0.235,
            0.395,
            [peak("low_q_broad", "OOP broad low-q component", 0.306, 0.285, 0.340, 0.070, 0.045, 0.150, 11, temperatures)],
        ),
        window(
            "feature_0430_early",
            "OOP early feature near q=0.430",
            0.395,
            0.470,
            [peak("feature_0430_early", "OOP early component near q=0.430", 0.430, 0.416, 0.442, 0.018, 0.006, 0.035, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_0530",
            "OOP feature near q=0.530",
            0.480,
            0.590,
            [peak("feature_0530", "OOP component near q=0.530", 0.530, 0.515, 0.545, 0.032, 0.012, 0.085, 11, temperatures)],
            temperature_min_C=50.0,
        ),
        window(
            "feature_0741_early",
            "OOP early feature near q=0.741",
            0.700,
            0.790,
            [peak("feature_0741_early", "OOP early component near q=0.741", 0.741, 0.732, 0.753, 0.020, 0.006, 0.040, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_1482_early",
            "OOP early feature near q=1.482",
            1.440,
            1.540,
            [peak("feature_1482_early", "OOP early component near q=1.482", 1.482, 1.468, 1.500, 0.018, 0.006, 0.045, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_1730_broad_early",
            "OOP provisional early broad shoulder near q=1.73",
            1.620,
            1.840,
            [peak("feature_1730_broad_early", "OOP early broad shoulder near q=1.73", 1.730, 1.670, 1.790, 0.105, 0.040, 0.220, 0, temperatures)],
            temperature_max_C=60.0,
        ),
        window(
            "feature_2217_early",
            "OOP early feature near q=2.217",
            2.160,
            2.250,
            [peak("feature_2217_early", "OOP early component near q=2.217", 2.217, 2.195, 2.225, 0.030, 0.006, 0.050, 0, temperatures)],
            temperature_max_C=60.0,
        ),
    ]


def validate_config(config: dict[str, Any]) -> None:
    windows = config["windows"]
    ordered = sorted(windows, key=lambda item: item["q_min"])
    for previous, current in zip(ordered, ordered[1:]):
        if previous["q_max"] > current["q_min"]:
            raise ValueError(
                f"Overlapping windows: {previous['key']} and {current['key']}"
            )
    temperatures = config["frame_temperatures_C"]
    for item in windows:
        if not item["q_min"] < item["q_max"]:
            raise ValueError(f"Invalid window bounds: {item['key']}")
        centre_ranges: list[tuple[float, float, str]] = []
        for component in item["peaks"]:
            if not (
                item["q_min"] < component["q_min"]
                < component["q_guess"]
                < component["q_max"] < item["q_max"]
            ):
                raise ValueError(f"Invalid q bounds: {component['key']}")
            if not (
                0 < component["fwhm_min"]
                <= component["fwhm_guess"]
                <= component["fwhm_max"]
            ):
                raise ValueError(f"Invalid FWHM bounds: {component['key']}")
            anchor = str(component["anchor_frame"])
            if anchor not in temperatures:
                raise ValueError(f"Missing anchor frame: {component['key']}")
            anchor_temperature = temperatures[anchor]
            if anchor_temperature < item.get("temperature_min_C", float("-inf")):
                raise ValueError(f"Anchor below temperature gate: {component['key']}")
            if anchor_temperature > item.get("temperature_max_C", float("inf")):
                raise ValueError(f"Anchor above temperature gate: {component['key']}")
            centre_ranges.append((component["q_min"], component["q_max"], component["key"]))
        centre_ranges.sort()
        for previous, current in zip(centre_ranges, centre_ranges[1:]):
            if previous[1] >= current[0]:
                raise ValueError(
                    f"Overlapping component centre bounds in {item['key']}: "
                    f"{previous[2]} and {current[2]}"
                )


def main() -> None:
    manifest = pd.read_csv(MANIFEST_PATH)
    manifest["frame"] = pd.to_numeric(manifest["frame"], errors="raise").astype(int)
    selected = manifest.loc[manifest["key_frame_selected"].astype(bool), "frame"].tolist()
    if selected != KEY_FRAMES:
        raise ValueError(f"Confirmed key frames changed: {selected}")
    temperatures = {
        str(int(row.frame)): float(row.measured_temperature_C)
        for row in manifest.itertuples(index=False)
    }
    builders = {"FR": fr_windows, "IP": ip_windows, "OOP": oop_windows}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for cut, build_windows in builders.items():
        config = {
            "frame_temperatures_C": temperatures,
            "profiles_to_compare": PROFILES,
            "fit_settings": FIT_SETTINGS,
            "windows": build_windows(temperatures),
            "analysis_scope": (
                "Scan-specific candidate-v2 key-frame screening for peak "
                "presence, q0/d and integrated area. No CCL interpretation."
            ),
            "representative_frames": KEY_FRAMES,
            "series_metadata": {
                "scan": int(SCAN),
                "sample_label": SAMPLE_LABEL,
                "cut": cut,
                "series_type": "secondary_in_situ_time_series",
                "frame_manifest": str(MANIFEST_PATH.resolve()),
                "key_frames": KEY_FRAMES,
                "visual_review": str(VISUAL_REVIEW_PATH.resolve()),
                "config_version": "candidate_v2",
                "late_cooling_frame": 35,
                "late_cooling_temperature_C": temperatures["35"],
                "anchor_policy": (
                    "Persistent/high-temperature features use frame 11; "
                    "features confined to frames 0-1 use frame 0."
                ),
                "scientific_status": (
                    "Windows are scan-specific visual candidates. They require "
                    "checkpoint overlay/residual review before selected-frame fitting."
                ),
            },
        }
        validate_config(config)
        output = OUTPUT_DIR / f"scan{SCAN}_{cut}_candidate_keyframes.json"
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
