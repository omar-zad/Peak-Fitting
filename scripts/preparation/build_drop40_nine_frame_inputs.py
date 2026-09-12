#!/usr/bin/env python3
"""Build reproducible nine-frame FR/IP/OOP fitting inputs for Drop40."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]


SCAN_TEMPERATURES_NOMINAL_C = {
    587193: 40.0,
    587194: 50.0,
    587195: 55.0,
    587196: 60.0,
    587197: 65.0,
    587198: 70.0,
    587199: 75.0,
    587200: 80.0,
    587201: 150.0,
}

CUTS = {
    "FR": "FULL",
    "IP": "IP",
    "OOP": "OOP",
}


def parse_dat_metadata(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    temperature_matches = re.findall(r"^temp3=([0-9.eE+-]+)", text, re.MULTILINE)
    if not temperature_matches:
        raise ValueError(f"No temp3 entry found in {path}")

    lines = text.splitlines()
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if {"transmission", "count_time", "norm"}.issubset(
                set(line.split("\t"))
            )
        ),
        None,
    )
    if header_index is None:
        raise ValueError(f"No named scan-data header found in {path}")
    columns = lines[header_index].split("\t")
    table_lines = [
        line
        for line in lines[header_index + 1 :]
        if re.match(r"^[0-9.+-]", line) and "\t" in line
    ]
    if not table_lines:
        raise ValueError(f"No scan-data row found after the header in {path}")
    fields = table_lines[-1].split("\t")
    if len(fields) != len(columns):
        raise ValueError(
            f"Scan-data header/row length mismatch in {path}: "
            f"{len(columns)} columns versus {len(fields)} values"
        )
    row = dict(zip(columns, fields))
    metadata = {
        "measured_temperature_C": float(temperature_matches[-1]),
        "transmission": float(row["transmission"]),
        "count_time_s": float(row["count_time"]),
        "norm_counts": float(row["norm"]),
    }
    if not all(np.isfinite(value) for value in metadata.values()):
        raise ValueError(f"Non-finite metadata value found in {path}")
    if metadata["norm_counts"] <= 0:
        raise ValueError(f"norm must be positive in {path}")
    if metadata["count_time_s"] <= 0:
        raise ValueError(f"count_time must be positive in {path}")
    return metadata


def locate_dat(data_root: Path, scan: int) -> Path:
    matches = sorted((data_root / "Data").glob(f"**/{scan}.dat"))
    if not matches:
        raise FileNotFoundError(f"Could not locate metadata file for scan {scan}")
    return matches[0]


def linecut_path(data_root: Path, scan: int, source_cut: str) -> Path:
    return (
        data_root
        / "Images"
        / f"scan_{scan}"
        / "Individual_1D"
        / "IP_OOP_manual_sector"
        / f"pilatus2-{scan}_{source_cut}_manual_wedge_1Dintegrations_Full.csv"
    )


def build_inputs(data_root: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    scans = list(SCAN_TEMPERATURES_NOMINAL_C)
    manifest_rows: list[dict[str, object]] = []

    metadata_by_scan: dict[int, dict[str, float]] = {}
    for scan in scans:
        dat_path = locate_dat(data_root, scan)
        metadata = parse_dat_metadata(dat_path)
        metadata_by_scan[scan] = metadata
        manifest_rows.append(
            {
                "frame": scan - 587000,
                "scan": scan,
                "nominal_temperature_C": SCAN_TEMPERATURES_NOMINAL_C[scan],
                **metadata,
                "metadata_file": str(dat_path),
            }
        )

    median_norm = float(
        np.median(
            [metadata_by_scan[scan]["norm_counts"] for scan in scans]
        )
    )
    for row in manifest_rows:
        row["normalization_reference_median_norm_counts"] = median_norm
        row["intensity_scale_median_norm_over_scan_norm"] = (
            median_norm / float(row["norm_counts"])
        )

    reference_q: np.ndarray | None = None
    for output_cut, source_cut in CUTS.items():
        profiles: dict[int, np.ndarray] = {}
        sources: dict[int, str] = {}
        local_q: np.ndarray | None = None

        for scan in scans:
            path = linecut_path(data_root, scan, source_cut)
            if not path.exists():
                raise FileNotFoundError(path)
            table = pd.read_csv(path)
            if table.shape[1] < 2:
                raise ValueError(f"Expected q and intensity columns in {path}")
            q = pd.to_numeric(table.iloc[:, 0], errors="coerce").to_numpy(float)
            intensity = pd.to_numeric(
                table.iloc[:, 1], errors="coerce"
            ).to_numpy(float)

            if local_q is None:
                local_q = q
            elif not np.array_equal(local_q, q, equal_nan=True):
                raise ValueError(
                    f"q grid for scan {scan} {output_cut} does not match "
                    "the first scan"
                )
            profiles[scan] = intensity
            sources[scan] = str(path)

        assert local_q is not None
        if reference_q is None:
            reference_q = local_q
        elif not np.array_equal(reference_q, local_q, equal_nan=True):
            raise ValueError(f"{output_cut} q grid differs from the FR q grid")

        output = pd.DataFrame({"Q A-1": local_q})
        for scan in scans:
            output[str(scan - 587000)] = profiles[scan]

        finite = np.isfinite(output["Q A-1"].to_numpy(float))
        for scan in scans:
            finite &= np.isfinite(output[str(scan - 587000)].to_numpy(float))
        output = output.loc[finite].reset_index(drop=True)

        output_path = output_dir / f"drop40_{output_cut}_9frames.txt"
        output.to_csv(
            output_path,
            sep="\t",
            index=False,
            float_format="%.10g",
        )

        for row in manifest_rows:
            scan = int(row["scan"])
            row[f"{output_cut}_source_file"] = sources[scan]
        print(
            f"{output_cut}: wrote {output_path} "
            f"({len(output)} q points, {len(scans)} frames)"
        )

        normalized = output.copy()
        for scan in scans:
            column = str(scan - 587000)
            scale = median_norm / metadata_by_scan[scan]["norm_counts"]
            normalized[column] = normalized[column] * scale
        normalized_path = (
            output_dir / f"drop40_{output_cut}_9frames_norm.txt"
        )
        normalized.to_csv(
            normalized_path,
            sep="\t",
            index=False,
            float_format="%.10g",
        )
        print(
            f"{output_cut}: wrote {normalized_path} "
            f"(scaled by median(norm)/scan norm)"
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "drop40_9frame_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"Manifest: wrote {manifest_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("external_data"),
        help="Root containing Data/ and Images/scan_*/ folders",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "local_results" / "drop40_inputs",
        help="Destination for the three input tables and manifest",
    )
    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    build_inputs(arguments.data_root.resolve(), arguments.output_dir.resolve())
