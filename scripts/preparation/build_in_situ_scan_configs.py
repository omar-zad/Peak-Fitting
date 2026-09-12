#!/usr/bin/env python3
"""Build scan-specific peak-fit configs from a frame manifest.

The production peak windows remain unchanged so that samples are compared on
the same scientific basis. Frame temperatures and anchor temperatures are
adapted to the measured in-situ series.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_CONFIGS = {
    "FR": PROJECT_ROOT / "configs" / "drop40" / "peakfit_config_FR_revised.json",
    "IP": PROJECT_ROOT / "configs" / "drop40" / "peakfit_config_IP_revised.json",
    "OOP": PROJECT_ROOT / "configs" / "drop40" / "peakfit_config_OOP_revised.json",
}


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"frame", "measured_temperature_C", "phase"}
    if not rows or not required.issubset(rows[0]):
        missing = sorted(required - set(rows[0] if rows else []))
        raise ValueError(f"Manifest is empty or missing columns: {', '.join(missing)}")
    frames = [int(row["frame"]) for row in rows]
    if frames != list(range(len(rows))):
        raise ValueError("Manifest frames must be continuous and zero-based.")
    return rows


def build_config(
    source_path: Path,
    manifest_path: Path,
    rows: list[dict[str, str]],
    scan: int,
    cut: str,
) -> dict:
    config = json.loads(source_path.read_text(encoding="utf-8"))
    temperatures = {
        row["frame"]: float(row["measured_temperature_C"]) for row in rows
    }
    low_anchor_frame = min(temperatures, key=temperatures.get)
    high_anchor_frame = max(temperatures, key=temperatures.get)
    low_anchor_temperature = temperatures[low_anchor_frame]
    high_anchor_temperature = temperatures[high_anchor_frame]

    config["frame_temperatures_C"] = temperatures
    for window in config["windows"]:
        for peak in window["peaks"]:
            old_anchor = float(peak["anchor_temperature_C"])
            peak["anchor_temperature_C"] = (
                low_anchor_temperature
                if old_anchor <= 60
                else high_anchor_temperature
            )

    config["series_metadata"] = {
        "scan": scan,
        "cut": cut,
        "series_type": "in_situ_time_series",
        "frame_manifest": str(manifest_path.resolve()),
        "frame_order": "zero-based acquisition order",
        "low_temperature_anchor": {
            "frame": int(low_anchor_frame),
            "measured_temperature_C": low_anchor_temperature,
        },
        "high_temperature_anchor": {
            "frame": int(high_anchor_frame),
            "measured_temperature_C": high_anchor_temperature,
        },
        "source_production_config": str(source_path.resolve()),
        "note": (
            "Peak windows and bounds are inherited unchanged for cross-sample "
            "comparability; only the frame-temperature map and anchor "
            "temperatures are adapted to this measured series."
        ),
    }
    return config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", type=int, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    rows = load_manifest(manifest_path)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for cut, source_path in SOURCE_CONFIGS.items():
        config = build_config(
            source_path=source_path,
            manifest_path=manifest_path,
            rows=rows,
            scan=args.scan,
            cut=cut,
        )
        output_path = output_dir / f"scan{args.scan}_{cut}.json"
        output_path.write_text(
            json.dumps(config, indent=2) + "\n",
            encoding="utf-8",
        )
        print(output_path)


if __name__ == "__main__":
    main()
