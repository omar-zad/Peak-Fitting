#!/usr/bin/env python3
"""Build visually revised key-frame candidate configs for scan 587250."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "configs" / "secondary_samples" / "scan_587250" / "candidate_v1"
OUTPUT = ROOT / "configs" / "secondary_samples" / "scan_587250" / "candidate_v2"


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


WINDOWS = {
    "FR": [
        window(
            "low_q_pair",
            "FR low-q broad maximum and higher-q shoulder",
            0.235,
            0.365,
            [
                peak("low_q_broad", "FR broad low-q component", 0.306, 0.282, 0.318, 0.060, 0.025, 0.120, 75, 72.3),
                peak("low_q_sharp", "FR sharper low-q shoulder", 0.326, 0.319, 0.338, 0.014, 0.008, 0.035, 30, 55.4),
            ],
        ),
        window(
            "feature_0430_early",
            "FR early feature near q = 0.43 A^-1",
            0.39,
            0.47,
            [peak("feature_0430_early", "FR early component near 0.43", 0.428, 0.414, 0.442, 0.017, 0.008, 0.040, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_0520",
            "FR evolving feature near q = 0.52 A^-1",
            0.47,
            0.575,
            [peak("feature_0520", "FR component near 0.52", 0.520, 0.505, 0.550, 0.035, 0.012, 0.090, 75, 72.3)],
        ),
        window(
            "feature_0738_early",
            "FR early feature near q = 0.74 A^-1",
            0.69,
            0.785,
            [peak("feature_0738_early", "FR early component near 0.74", 0.738, 0.722, 0.752, 0.017, 0.008, 0.045, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_1500_1535_early",
            "FR early q = 1.50 / 1.53 A^-1 region",
            1.45,
            1.575,
            [
                peak("feature_1500_early", "FR early lower-q component", 1.500, 1.475, 1.514, 0.030, 0.008, 0.070, 30, 55.4),
                peak("feature_1535_early", "FR early higher-q component", 1.535, 1.515, 1.552, 0.018, 0.008, 0.052, 0, 32.3),
            ],
            temperature_max=65.0,
        ),
        window(
            "feature_1535_late",
            "FR intermediate/high-temperature feature near q = 1.535 A^-1",
            1.49,
            1.585,
            [peak("feature_1535_late", "FR later component near 1.535", 1.535, 1.518, 1.552, 0.022, 0.008, 0.060, 75, 72.3)],
            temperature_min=65.0,
        ),
        window(
            "feature_1680_early",
            "FR early feature near q = 1.68 A^-1",
            1.62,
            1.735,
            [peak("feature_1680_early", "FR early component near 1.68", 1.680, 1.650, 1.700, 0.022, 0.008, 0.060, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_1720_late",
            "FR later broad feature near q = 1.72 A^-1",
            1.67,
            1.79,
            [peak("feature_1720_late", "FR later broad component near 1.72", 1.720, 1.700, 1.745, 0.045, 0.014, 0.110, 75, 72.3)],
            temperature_min=65.0,
        ),
        window(
            "feature_1960_early",
            "FR early feature near q = 1.96 A^-1",
            1.90,
            2.015,
            [peak("feature_1960_early", "FR early component near 1.96", 1.960, 1.930, 1.990, 0.025, 0.008, 0.065, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_2067_early",
            "FR early feature near q = 2.067 A^-1",
            2.02,
            2.105,
            [peak("feature_2067_early", "FR early component near 2.067", 2.067, 2.050, 2.085, 0.020, 0.008, 0.050, 30, 55.4)],
            temperature_max=65.0,
        ),
        window(
            "feature_2155",
            "FR strong feature near q = 2.155 A^-1",
            2.08,
            2.22,
            [peak("feature_2155", "FR strong component near 2.155", 2.155, 2.135, 2.172, 0.040, 0.012, 0.095, 75, 72.3)],
        ),
    ],
    "IP": [
        window(
            "low_q_pair",
            "IP low-q broad maximum and q = 0.324 A^-1 component",
            0.235,
            0.365,
            [
                peak("low_q_broad", "IP broad low-q component", 0.292, 0.265, 0.315, 0.060, 0.025, 0.125, 75, 72.3),
                peak("low_q_sharp", "IP sharper component near 0.324", 0.324, 0.316, 0.337, 0.015, 0.008, 0.038, 75, 72.3),
            ],
        ),
        window(
            "feature_0415_early",
            "IP early feature near q = 0.415 A^-1",
            0.375,
            0.45,
            [peak("feature_0415_early", "IP early component near 0.415", 0.415, 0.398, 0.430, 0.018, 0.008, 0.045, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_0448_0502_late",
            "IP later q = 0.448 / 0.502 A^-1 pair",
            0.415,
            0.545,
            [
                peak("feature_0448_late", "IP later lower-q component", 0.448, 0.435, 0.462, 0.025, 0.008, 0.060, 75, 72.3),
                peak("feature_0502_late", "IP later higher-q component", 0.502, 0.486, 0.518, 0.030, 0.010, 0.070, 129, 83.3),
            ],
            temperature_min=60.0,
        ),
        window(
            "feature_0730_early",
            "IP early feature near q = 0.73 A^-1",
            0.685,
            0.78,
            [peak("feature_0730_early", "IP early component near 0.73", 0.730, 0.715, 0.748, 0.018, 0.008, 0.045, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_1505_1538_early",
            "IP early q = 1.50 / 1.54 A^-1 region",
            1.46,
            1.575,
            [
                peak("feature_1505_early", "IP early lower-q component", 1.505, 1.482, 1.518, 0.030, 0.008, 0.070, 30, 55.4),
                peak("feature_1538_early", "IP early higher-q component", 1.538, 1.519, 1.557, 0.020, 0.008, 0.055, 11, 30.0),
            ],
            temperature_max=65.0,
        ),
        window(
            "feature_1538_late",
            "IP later feature near q = 1.54 A^-1",
            1.49,
            1.585,
            [peak("feature_1538_late", "IP later component near 1.54", 1.540, 1.522, 1.558, 0.025, 0.008, 0.065, 129, 83.3)],
            temperature_min=65.0,
        ),
        window(
            "feature_1680_early",
            "IP early feature near q = 1.68 A^-1",
            1.62,
            1.735,
            [peak("feature_1680_early", "IP early component near 1.68", 1.680, 1.650, 1.705, 0.024, 0.008, 0.065, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_1728_late",
            "IP later broad feature near q = 1.73 A^-1",
            1.67,
            1.80,
            [peak("feature_1728_late", "IP later broad component near 1.73", 1.728, 1.705, 1.752, 0.050, 0.014, 0.120, 75, 72.3)],
            temperature_min=65.0,
        ),
        window(
            "feature_1980_early",
            "IP early feature near q = 1.98 A^-1",
            1.92,
            2.025,
            [peak("feature_1980_early", "IP early component near 1.98", 1.980, 1.955, 2.000, 0.024, 0.008, 0.060, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_2085_early",
            "IP early feature near q = 2.085 A^-1",
            2.035,
            2.125,
            [peak("feature_2085_early", "IP early component near 2.085", 2.085, 2.065, 2.102, 0.022, 0.008, 0.055, 30, 55.4)],
            temperature_max=65.0,
        ),
        window(
            "feature_2186",
            "IP strong feature near q = 2.186 A^-1",
            2.11,
            2.245,
            [peak("feature_2186", "IP strong component near 2.186", 2.186, 2.168, 2.202, 0.040, 0.012, 0.095, 75, 72.3)],
        ),
    ],
    "OOP": [
        window(
            "low_q_broad",
            "OOP broad low-q maximum",
            0.235,
            0.385,
            [peak("low_q_broad", "OOP broad low-q component", 0.306, 0.282, 0.335, 0.075, 0.035, 0.145, 75, 72.3)],
        ),
        window(
            "feature_0430_early",
            "OOP early feature near q = 0.43 A^-1",
            0.39,
            0.47,
            [peak("feature_0430_early", "OOP early component near 0.43", 0.428, 0.414, 0.443, 0.017, 0.008, 0.045, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_0525",
            "OOP evolving feature near q = 0.525 A^-1",
            0.47,
            0.59,
            [peak("feature_0525", "OOP component near 0.525", 0.525, 0.505, 0.555, 0.035, 0.012, 0.090, 75, 72.3)],
        ),
        window(
            "feature_0738_early",
            "OOP early feature near q = 0.74 A^-1",
            0.69,
            0.785,
            [peak("feature_0738_early", "OOP early component near 0.74", 0.738, 0.720, 0.752, 0.017, 0.008, 0.045, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_1465_1532_early",
            "OOP early q = 1.465 / 1.532 A^-1 region",
            1.42,
            1.575,
            [
                peak("feature_1465_early", "OOP early lower-q component", 1.465, 1.445, 1.485, 0.030, 0.008, 0.070, 22, 29.6),
                peak("feature_1532_early", "OOP early higher-q component", 1.532, 1.512, 1.550, 0.022, 0.008, 0.060, 0, 32.3),
            ],
            temperature_max=65.0,
        ),
        window(
            "feature_1532_late",
            "OOP later feature near q = 1.532 A^-1",
            1.48,
            1.585,
            [peak("feature_1532_late", "OOP later component near 1.532", 1.532, 1.515, 1.550, 0.025, 0.008, 0.065, 75, 72.3)],
            temperature_min=65.0,
        ),
        window(
            "feature_1680_early",
            "OOP early feature near q = 1.68 A^-1",
            1.62,
            1.735,
            [peak("feature_1680_early", "OOP early component near 1.68", 1.680, 1.650, 1.705, 0.022, 0.008, 0.065, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_1720_late",
            "OOP later broad feature near q = 1.72 A^-1",
            1.67,
            1.79,
            [peak("feature_1720_late", "OOP later broad component near 1.72", 1.720, 1.700, 1.745, 0.045, 0.014, 0.110, 75, 72.3)],
            temperature_min=65.0,
        ),
        window(
            "feature_1950_early",
            "OOP early feature near q = 1.95 A^-1",
            1.89,
            2.015,
            [peak("feature_1950_early", "OOP early component near 1.95", 1.950, 1.925, 1.980, 0.025, 0.008, 0.065, 22, 29.6)],
            temperature_max=65.0,
        ),
        window(
            "feature_2055_early",
            "OOP early feature near q = 2.055 A^-1",
            2.01,
            2.105,
            [peak("feature_2055_early", "OOP early component near 2.055", 2.055, 2.035, 2.075, 0.022, 0.008, 0.055, 30, 55.4)],
            temperature_max=65.0,
        ),
        window(
            "feature_2145",
            "OOP strong feature near q = 2.145 A^-1",
            2.07,
            2.215,
            [peak("feature_2145", "OOP strong component near 2.145", 2.145, 2.125, 2.162, 0.042, 0.012, 0.100, 75, 72.3)],
        ),
    ],
}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for cut, windows in WINDOWS.items():
        source = SOURCE / f"scan587250_{cut}_candidate_keyframes.json"
        config = json.loads(source.read_text(encoding="utf-8"))
        config["windows"] = windows
        config["series_metadata"]["template_config"] = str(source.resolve())
        config["series_metadata"]["scientific_status"] = (
            "Candidate-v2 windows revised from scan 587250's own unsmoothed "
            "key-frame overlays; retain only peaks passing numerical and visual QC."
        )
        output = OUTPUT / f"scan587250_{cut}_candidate_keyframes.json"
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
