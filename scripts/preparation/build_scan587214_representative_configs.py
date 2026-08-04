#!/usr/bin/env python3
"""Build the seven-frame model-development configs for scan 587214."""

from __future__ import annotations

import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = (
    PROJECT_ROOT
    / "configs"
    / "in_situ"
    / "scan_587214"
    / "revised_candidate"
)
OUTPUT_DIR = (
    PROJECT_ROOT
    / "configs"
    / "in_situ"
    / "scan_587214"
    / "representative"
)
REPRESENTATIVE_FRAMES = [0, 25, 50, 57, 71, 85, 99]


def ip_low_q_window() -> dict:
    return {
        "key": "low_q_pair",
        "title": "IP unassigned low-q broad and sharp components",
        "q_min": 0.245,
        "q_max": 0.400,
        "peaks": [
            {
                "key": "low_q_broad",
                "label": "IP unassigned broad component near 0.30",
                "q_guess": 0.305,
                "q_min": 0.282,
                "q_max": 0.330,
                "fwhm_guess": 0.055,
                "fwhm_min": 0.025,
                "fwhm_max": 0.115,
                "anchor_frame": 85,
                "anchor_temperature_C": 76.3,
                "provisional_assignment": True,
            },
            {
                "key": "low_q_sharp",
                "label": "IP unassigned sharp component near 0.35",
                "q_guess": 0.350,
                "q_min": 0.338,
                "q_max": 0.362,
                "fwhm_guess": 0.020,
                "fwhm_min": 0.009,
                "fwhm_max": 0.055,
                "anchor_frame": 85,
                "anchor_temperature_C": 76.3,
                "provisional_assignment": True,
            },
        ],
    }


def build_cut(cut: str) -> Path:
    source = SOURCE_DIR / f"scan587214_{cut}_frame057_candidate.json"
    config = json.loads(source.read_text(encoding="utf-8"))
    config["analysis_scope"] = (
        "Representative-frame peak detection, q0, and integrated area; "
        "no CCL interpretation."
    )
    config["representative_frames"] = REPRESENTATIVE_FRAMES
    settings = config["fit_settings"]
    settings["multistarts_model_selection"] = 4
    settings["multistarts_final"] = 6

    if cut == "IP":
        config["windows"][0] = ip_low_q_window()

    for window in config["windows"]:
        window["title"] = window["title"].replace(
            "provisional", "unassigned"
        )
        for peak in window["peaks"]:
            peak["label"] = peak["label"].replace(
                "provisional", "unassigned"
            )
            peak.setdefault("anchor_frame", 57)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / f"scan587214_{cut}_keyframes.json"
    output.write_text(
        json.dumps(config, indent=2) + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    for cut in ("FR", "IP", "OOP"):
        print(build_cut(cut))


if __name__ == "__main__":
    main()
