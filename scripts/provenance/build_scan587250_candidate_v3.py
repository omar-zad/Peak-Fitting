#!/usr/bin/env python3
"""Refine scan 587250 candidate configs after the two-regime fit gate."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "configs" / "secondary_samples" / "scan_587250" / "candidate_v2"
OUTPUT = ROOT / "configs" / "secondary_samples" / "scan_587250" / "candidate_v3"


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
    anchor_temperature: float,
) -> dict:
    return {
        "key": key,
        "label": label,
        "q_guess": q_guess,
        "q_min": q_min,
        "q_max": q_max,
        "fwhm_guess": fwhm_guess,
        "fwhm_min": fwhm_min,
        "fwhm_max": fwhm_max,
        "anchor_temperature_C": anchor_temperature,
        "anchor_frame": anchor_frame,
        "provisional_assignment": True,
    }


def window(
    key: str,
    title: str,
    q_min: float,
    q_max: float,
    peaks: list[dict],
    *,
    temperature_min: float | None = None,
    temperature_max: float | None = None,
) -> dict:
    value = {
        "key": key,
        "title": title,
        "q_min": q_min,
        "q_max": q_max,
        "peaks": peaks,
        "candidate_window_requires_visual_review": True,
    }
    if temperature_min is not None:
        value["temperature_min_C"] = temperature_min
    if temperature_max is not None:
        value["temperature_max_C"] = temperature_max
    return value


FR_EARLY_MIDQ = window(
    "feature_0540_0580_early",
    "FR early q = 0.54 / 0.58 A^-1 structure",
    0.485,
    0.61,
    [
        peak("feature_0540_early", "FR early lower-q component", 0.540, 0.522, 0.560, 0.028, 0.008, 0.070, 22, 29.6),
        peak("feature_0580_early", "FR early higher-q component", 0.580, 0.565, 0.600, 0.022, 0.008, 0.060, 22, 29.6),
    ],
    temperature_max=65.0,
)

FR_LATE_MIDQ = window(
    "feature_0520_late",
    "FR later feature near q = 0.52 A^-1",
    0.47,
    0.575,
    [peak("feature_0520_late", "FR later component near 0.52", 0.520, 0.505, 0.535, 0.035, 0.012, 0.085, 75, 72.3)],
    temperature_min=65.0,
)

OOP_EARLY_MIDQ = window(
    "feature_0545_0585_early",
    "OOP early q = 0.545 / 0.585 A^-1 structure",
    0.49,
    0.615,
    [
        peak("feature_0545_early", "OOP early lower-q component", 0.545, 0.525, 0.563, 0.026, 0.008, 0.065, 22, 29.6),
        peak("feature_0585_early", "OOP early higher-q component", 0.585, 0.570, 0.602, 0.020, 0.008, 0.055, 22, 29.6),
    ],
    temperature_max=65.0,
)

OOP_LATE_MIDQ = window(
    "feature_0525_late",
    "OOP later feature near q = 0.525 A^-1",
    0.47,
    0.59,
    [peak("feature_0525_late", "OOP later component near 0.525", 0.525, 0.505, 0.540, 0.035, 0.012, 0.085, 75, 72.3)],
    temperature_min=65.0,
)

IP_EARLY_HIGHQ_PAIR = window(
    "feature_1655_1682_early",
    "IP early q = 1.655 / 1.682 A^-1 structure",
    1.61,
    1.735,
    [
        peak("feature_1655_early", "IP early lower-q shoulder", 1.655, 1.638, 1.668, 0.020, 0.008, 0.050, 22, 29.6),
        peak("feature_1682_early", "IP early higher-q component", 1.682, 1.669, 1.705, 0.024, 0.008, 0.060, 22, 29.6),
    ],
    temperature_max=65.0,
)


def replace_window(windows: list[dict], key: str, replacements: list[dict]) -> list[dict]:
    revised: list[dict] = []
    for item in windows:
        if item["key"] == key:
            revised.extend(replacements)
        else:
            revised.append(item)
    return revised


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for cut in ("FR", "IP", "OOP"):
        source = SOURCE / f"scan587250_{cut}_candidate_keyframes.json"
        config = json.loads(source.read_text(encoding="utf-8"))
        if cut == "FR":
            config["windows"] = replace_window(
                config["windows"], "feature_0520", [FR_EARLY_MIDQ, FR_LATE_MIDQ]
            )
        elif cut == "IP":
            config["windows"] = replace_window(
                config["windows"], "feature_1680_early", [IP_EARLY_HIGHQ_PAIR]
            )
        else:
            config["windows"] = replace_window(
                config["windows"], "feature_0525", [OOP_EARLY_MIDQ, OOP_LATE_MIDQ]
            )
        for item in config["windows"]:
            if item["key"] in {"feature_2155", "feature_2186", "feature_2145"}:
                item["reporting_scope"] = "q_position_only_reference_or_substrate_candidate"
                item["interpretation_note"] = (
                    "Strong high-q feature retained for position tracking, but its "
                    "asymmetric line shape and possible substrate origin require manual QC."
                )
        config["series_metadata"]["template_config"] = str(source.resolve())
        config["series_metadata"]["scientific_status"] = (
            "Candidate-v3 incorporates the low/high-temperature fit-gate review. "
            "Strong asymmetric high-q reference/substrate candidates require manual QC."
        )
        output = OUTPUT / f"scan587250_{cut}_candidate_keyframes.json"
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
