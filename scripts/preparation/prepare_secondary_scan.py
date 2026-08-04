#!/usr/bin/env python3
"""Prepare one secondary-sample GIWAXS scan for staged peak fitting.

The script validates the temperature table, beamline monitor metadata, and
FR/IP/OOP line cuts before writing reproducible raw and d5i-normalized fitting
inputs. It also proposes a small phase-aware set of key frames. No fitting is
performed.
"""

from __future__ import annotations

import argparse
import io
import json
import posixpath
import re
import zipfile
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REGISTRY = (
    PROJECT_ROOT / "configs" / "secondary_samples" / "scan_registry.json"
)
CUTS = ("FR", "IP", "OOP")


def parse_frame_list(value: str | None) -> list[int] | None:
    if value is None:
        return None
    frames = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not frames:
        raise argparse.ArgumentTypeError("Provide at least one frame number.")
    if len(frames) != len(set(frames)):
        raise argparse.ArgumentTypeError("Key-frame list contains duplicates.")
    return frames


def load_registry(path: Path, scan: str) -> dict[str, Any]:
    registry = json.loads(path.read_text(encoding="utf-8"))
    try:
        entry = registry["scans"][scan]
    except KeyError as exc:
        raise KeyError(f"Scan {scan} is absent from {path}.") from exc
    required = {
        "sample_label",
        "sample_slug",
        "temperature_file",
        "dat_file",
        "profile_files",
    }
    missing = sorted(required - set(entry))
    if missing:
        raise ValueError(
            f"Scan {scan} registry entry is missing: {', '.join(missing)}"
        )
    missing_cuts = sorted(set(CUTS) - set(entry["profile_files"]))
    if missing_cuts:
        raise ValueError(
            f"Scan {scan} is missing profile paths for: {', '.join(missing_cuts)}"
        )
    return entry


def _excel_column_index(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference.upper())
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {reference!r}")
    index = 0
    for character in match.group(1):
        index = index * 26 + (ord(character) - ord("A") + 1)
    return index - 1


def _read_simple_xlsx(path: Path, sheet_name: str | None = None) -> pd.DataFrame:
    """Read one worksheet without requiring an optional Excel engine."""
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    office_rel_ns = (
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    )
    package_rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"
    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{main_ns}}}si"):
                text = "".join(
                    node.text or "" for node in item.iter(f"{{{main_ns}}}t")
                )
                shared_strings.append(text)
        workbook_root = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet_nodes = workbook_root.findall(f".//{{{main_ns}}}sheet")
        if not sheet_nodes:
            raise ValueError(f"{path.name} contains no worksheets.")
        relationships_root = ET.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        relationships = {
            node.attrib["Id"]: node
            for node in relationships_root.findall(
                f"{{{package_rel_ns}}}Relationship"
            )
        }
        worksheet_nodes = []
        for node in sheet_nodes:
            relationship_id = node.attrib.get(f"{{{office_rel_ns}}}id")
            relationship = relationships.get(relationship_id)
            if relationship is None:
                continue
            relationship_type = relationship.attrib.get("Type", "")
            target = relationship.attrib.get("Target", "")
            if relationship_type.endswith("/worksheet") or "worksheets/" in target:
                worksheet_nodes.append((node, relationship))
        if not worksheet_nodes:
            raise ValueError(f"{path.name} contains no data worksheets.")
        available_sheets = [node.attrib["name"] for node, _ in worksheet_nodes]
        if sheet_name is None:
            selected_sheet, relationship = worksheet_nodes[0]
        else:
            selected = next(
                (
                    (node, related)
                    for node, related in worksheet_nodes
                    if node.attrib["name"] == sheet_name
                ),
                None,
            )
            if selected is None:
                raise ValueError(
                    f"Worksheet {sheet_name!r} is absent from {path.name}; "
                    f"available sheets: {available_sheets}"
                )
            selected_sheet, relationship = selected
        target = relationship.attrib["Target"]
        worksheet_path = (
            target.lstrip("/")
            if target.startswith("/")
            else posixpath.normpath(posixpath.join("xl", target))
        )
        root = ET.fromstring(archive.read(worksheet_path))
        rows: list[list[object]] = []
        for row in root.iter(f"{{{main_ns}}}row"):
            values: dict[int, object] = {}
            for cell in row.findall(f"{{{main_ns}}}c"):
                column = _excel_column_index(cell.attrib.get("r", "A1"))
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{main_ns}}}v")
                if cell_type == "inlineStr":
                    value: object = "".join(
                        node.text or ""
                        for node in cell.iter(f"{{{main_ns}}}t")
                    )
                elif value_node is None or value_node.text is None:
                    value = None
                elif cell_type == "s":
                    value = shared_strings[int(value_node.text)]
                else:
                    raw = value_node.text
                    try:
                        value = float(raw)
                    except ValueError:
                        value = raw
                values[column] = value
            if values:
                width = max(values) + 1
                rows.append([values.get(index) for index in range(width)])
    if not rows:
        raise ValueError(
            f"Worksheet {selected_sheet.attrib['name']!r} in {path.name} is empty."
        )
    width = max(len(row) for row in rows)
    padded = [row + [None] * (width - len(row)) for row in rows]
    headers = [
        str(int(value))
        if isinstance(value, float) and value.is_integer()
        else str(value)
        if value is not None
        else ""
        for value in padded[0]
    ]
    return pd.DataFrame(padded[1:], columns=headers)


def load_temperature_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        table = pd.read_csv(path)
    elif path.suffix.lower() in {".xlsx", ".xls"}:
        try:
            table = pd.read_excel(path)
        except ImportError:
            if path.suffix.lower() != ".xlsx":
                raise
            table = _read_simple_xlsx(path)
    else:
        raise ValueError(f"Unsupported temperature-table format: {path.suffix}")

    required = {"frame", "time_s", "temp3_C"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(
            f"{path.name} is missing columns: {', '.join(missing)}"
        )
    result = table.loc[:, ["frame", "time_s", "temp3_C"]].copy()
    result.columns = ["frame", "acquisition_time_s", "measured_temperature_C"]
    result["frame"] = pd.to_numeric(result["frame"], errors="raise").astype(int)
    for column in ("acquisition_time_s", "measured_temperature_C"):
        result[column] = pd.to_numeric(result[column], errors="raise").astype(float)
    expected = list(range(len(result)))
    if result["frame"].tolist() != expected:
        raise ValueError(
            "Temperature-table frames must be continuous, zero-based, and in "
            "acquisition order."
        )
    if not np.isfinite(result[["acquisition_time_s", "measured_temperature_C"]]).all().all():
        raise ValueError("Temperature table contains non-finite time/temperature values.")
    return result


def load_srs_table(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.startswith("testMotor1\t") and "frameNo" in line and "d5i" in line
        ),
        None,
    )
    if header_index is None:
        raise ValueError(f"Could not locate the SRS data table in {path}.")
    table = pd.read_csv(
        io.StringIO("\n".join(lines[header_index:])),
        sep="\t",
        engine="python",
    )
    table.columns = [str(column).strip() for column in table.columns]
    required = {"frameNo", "Temperature", "Timer", "d5i", "count_time", "transmission"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{path.name} SRS table is missing: {', '.join(missing)}")
    table["frameNo"] = pd.to_numeric(table["frameNo"], errors="coerce")
    table = table.loc[table["frameNo"].notna()].copy()
    for column in required:
        table[column] = pd.to_numeric(table[column], errors="raise")
    table["frame"] = table["frameNo"].astype(int) - 1
    if table["frame"].tolist() != list(range(len(table))):
        raise ValueError(
            f"{path.name} frameNo values do not map continuously to zero-based frames."
        )
    return table


def load_profile(
    path: Path,
    expected_frames: list[int],
    sheet_name: str | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    source_label = path.name + (f" [{sheet_name}]" if sheet_name else "")
    if path.suffix.lower() == ".csv":
        if sheet_name is not None:
            raise ValueError(f"A worksheet was specified for CSV source {path.name}.")
        table = pd.read_csv(path)
    elif path.suffix.lower() in {".xlsx", ".xls"}:
        try:
            table = pd.read_excel(path, sheet_name=sheet_name or 0)
        except ImportError:
            if path.suffix.lower() != ".xlsx":
                raise
            table = _read_simple_xlsx(path, sheet_name=sheet_name)
    else:
        raise ValueError(f"Unsupported profile format: {path.suffix}")
    if table.shape[1] != len(expected_frames) + 1:
        raise ValueError(
            f"{source_label} has {table.shape[1] - 1} frame columns; "
            f"expected {len(expected_frames)}."
        )
    q = pd.to_numeric(table.iloc[:, 0], errors="raise").to_numpy(float)
    if len(q) < 8 or not np.isfinite(q).all() or not np.all(np.diff(q) > 0):
        raise ValueError(f"{source_label} q grid must be finite and strictly increasing.")
    frame_columns = []
    for column in table.columns[1:]:
        try:
            frame_columns.append(int(str(column).strip()))
        except ValueError as exc:
            raise ValueError(
                f"Cannot interpret frame column {column!r} in {source_label}."
            ) from exc
    if frame_columns != expected_frames:
        raise ValueError(
            f"{source_label} frame columns do not match zero-based acquisition order."
        )
    intensity = table.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    intensity.columns = [str(frame) for frame in expected_frames]
    return q, intensity


def phase_boundaries(temperature: np.ndarray) -> tuple[int, int, dict[str, float]]:
    n = len(temperature)
    dynamic = float(np.nanmax(temperature) - np.nanmin(temperature))
    if n < 3 or dynamic < 2.0:
        return 0, n - 1, {
            "baseline_temperature_C": float(np.nanmedian(temperature)),
            "high_threshold_C": float(np.nanmax(temperature)),
            "dynamic_range_C": dynamic,
        }
    baseline_count = min(5, n)
    baseline = float(np.nanmedian(temperature[:baseline_count]))
    rise_threshold = baseline + max(1.0, 0.05 * dynamic)
    immediate_rise = float(temperature[1] - temperature[0]) >= max(
        1.0, 0.03 * dynamic
    )
    if immediate_rise:
        ramp_start = 0
        baseline = float(temperature[0])
        rise_threshold = baseline
    else:
        ramp_candidates = np.flatnonzero(temperature >= rise_threshold)
        ramp_start = int(ramp_candidates[0]) if ramp_candidates.size else 0
    high_threshold = float(np.nanmax(temperature) - max(1.0, 0.03 * dynamic))
    high_candidates = np.flatnonzero(
        (np.arange(n) >= ramp_start) & (temperature >= high_threshold)
    )
    high_start = int(high_candidates[0]) if high_candidates.size else int(np.nanargmax(temperature))
    high_start = max(high_start, ramp_start)
    return ramp_start, high_start, {
        "baseline_temperature_C": baseline,
        "rise_threshold_C": rise_threshold,
        "high_threshold_C": high_threshold,
        "dynamic_range_C": dynamic,
    }


def _add_key(
    chosen: dict[int, str], frame: int, role: str, frame_count: int
) -> None:
    frame = max(0, min(int(frame), frame_count - 1))
    chosen.setdefault(frame, role)


def suggest_key_frames(
    temperatures: np.ndarray,
    ramp_start: int,
    high_start: int,
) -> dict[int, str]:
    """Return at most seven phase-aware representative frames."""
    n = len(temperatures)
    chosen: dict[int, str] = {}
    _add_key(chosen, 0, "series_start", n)

    low_count = ramp_start
    if low_count >= 6:
        _add_key(chosen, (ramp_start - 1) // 2, "low_period_midpoint", n)
        _add_key(chosen, ramp_start - 1, "low_period_end", n)

    ramp_length = high_start - ramp_start + 1
    if ramp_length > 1:
        fractions = (0.25, 0.50, 0.75) if low_count < 6 else (0.50,)
        for fraction in fractions:
            frame = int(round(ramp_start + fraction * (ramp_length - 1)))
            _add_key(chosen, frame, f"heating_ramp_{int(fraction * 100)}pct", n)

    _add_key(chosen, high_start, "high_temperature_start", n)
    high_count = n - high_start
    if high_count >= 6:
        _add_key(chosen, high_start + (high_count - 1) // 2, "high_period_midpoint", n)
    _add_key(chosen, n - 1, "series_end", n)

    if len(chosen) > 7:
        priority = {
            "series_start": 0,
            "low_period_midpoint": 1,
            "low_period_end": 2,
            "heating_ramp_50pct": 3,
            "high_temperature_start": 4,
            "high_period_midpoint": 5,
            "series_end": 6,
        }
        retained = sorted(
            chosen.items(), key=lambda item: priority.get(item[1], 99)
        )[:7]
        chosen = dict(retained)
    return dict(sorted(chosen.items()))


def role_for_confirmed_frame(
    frame: int, ramp_start: int, high_start: int, frame_count: int
) -> str:
    if frame == 0:
        return "series_start"
    if frame == frame_count - 1:
        return "series_end"
    if frame < ramp_start:
        return "selected_low_temperature_period"
    if frame < high_start:
        return "selected_heating_ramp"
    if frame == high_start:
        return "selected_high_temperature_start"
    return "selected_high_temperature_period"


def prepare(args: argparse.Namespace) -> dict[str, Path]:
    registry_path = args.registry.expanduser().resolve()
    entry = load_registry(registry_path, args.scan)
    temperature_path = Path(entry["temperature_file"]).expanduser().resolve()
    dat_path = Path(entry["dat_file"]).expanduser().resolve()
    profile_sources: dict[str, dict[str, Any]] = {}
    for cut in CUTS:
        source = entry["profile_files"][cut]
        if isinstance(source, str):
            source_path = source
            source_sheet = None
        elif isinstance(source, dict) and "path" in source:
            source_path = source["path"]
            source_sheet = source.get("sheet")
        else:
            raise ValueError(
                f"Profile source for {cut} must be a path string or an object "
                "containing 'path' and optional 'sheet'."
            )
        profile_sources[cut] = {
            "path": Path(source_path).expanduser().resolve(),
            "sheet": source_sheet,
        }
    profile_paths = {cut: source["path"] for cut, source in profile_sources.items()}
    source_paths = [temperature_path, dat_path, *profile_paths.values()]
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing source files:\n  " + "\n  ".join(missing))

    temperature = load_temperature_table(temperature_path)
    srs = load_srs_table(dat_path)
    if len(temperature) != len(srs):
        raise ValueError(
            f"Temperature table has {len(temperature)} frames but {dat_path.name} "
            f"has {len(srs)}."
        )
    temperature_difference = np.abs(
        temperature["measured_temperature_C"].to_numpy(float)
        - srs["Temperature"].to_numpy(float)
    )
    if float(np.nanmax(temperature_difference)) > 0.25:
        raise ValueError(
            "Temperature table and SRS Temperature differ by more than 0.25 C; "
            "verify frame alignment before fitting."
        )

    frames = temperature["frame"].tolist()
    profiles: dict[str, tuple[np.ndarray, pd.DataFrame]] = {
        cut: load_profile(
            profile_sources[cut]["path"],
            frames,
            sheet_name=profile_sources[cut]["sheet"],
        )
        for cut in CUTS
    }
    q_reference = profiles["FR"][0]
    for cut in CUTS[1:]:
        q = profiles[cut][0]
        if q.shape != q_reference.shape or not np.allclose(
            q, q_reference, rtol=0.0, atol=1e-10, equal_nan=True
        ):
            raise ValueError(f"{cut} q grid does not match the FR q grid.")

    profile_coverage: dict[str, dict[str, float | int]] = {}
    for cut, (q, intensity) in profiles.items():
        finite = np.isfinite(intensity.to_numpy(float))
        finite_per_frame = finite.sum(axis=0)
        if np.any(finite_per_frame < 8):
            bad_frames = np.flatnonzero(finite_per_frame < 8).tolist()
            raise ValueError(
                f"{cut} has fewer than eight finite intensity points in frames "
                f"{bad_frames}."
            )
        first_finite = np.argmax(finite, axis=0)
        last_finite = len(q) - 1 - np.argmax(finite[::-1], axis=0)
        profile_coverage[cut] = {
            "nonfinite_intensity_cells": int(np.count_nonzero(~finite)),
            "minimum_finite_q_points_per_frame": int(finite_per_frame.min()),
            "all_frame_usable_q_min_Ainv": float(np.max(q[first_finite])),
            "all_frame_usable_q_max_Ainv": float(np.min(q[last_finite])),
        }

    monitor = srs["d5i"].to_numpy(float)
    if not np.isfinite(monitor).all() or np.any(monitor <= 0):
        raise ValueError("d5i monitor contains non-finite or non-positive values.")
    monitor_median = float(np.median(monitor))
    scale = monitor_median / monitor

    ramp_start, high_start, phase_details = phase_boundaries(
        temperature["measured_temperature_C"].to_numpy(float)
    )
    automatic_keys = suggest_key_frames(
        temperature["measured_temperature_C"].to_numpy(float),
        ramp_start,
        high_start,
    )
    registry_suggestions = entry.get("suggested_key_frames")
    if registry_suggestions is not None:
        invalid_suggestions = [
            int(frame) for frame in registry_suggestions if int(frame) not in frames
        ]
        if invalid_suggestions:
            raise ValueError(
                f"Registry key-frame suggestions are outside the scan: {invalid_suggestions}"
            )
        suggested_keys = {
            int(frame): role_for_confirmed_frame(
                int(frame), ramp_start, high_start, len(frames)
            )
            for frame in registry_suggestions
        }
        suggestion_source = "registry_phase_review"
    else:
        suggested_keys = automatic_keys
        suggestion_source = "automatic_temperature_trace"
    cli_keys = args.key_frames
    registry_keys = entry.get("key_frames")
    confirmed_frames = cli_keys if cli_keys is not None else registry_keys
    if confirmed_frames is not None:
        confirmed_frames = [int(frame) for frame in confirmed_frames]
        invalid = [frame for frame in confirmed_frames if frame not in frames]
        if invalid:
            raise ValueError(f"Confirmed key frames are outside the scan: {invalid}")
        selected_keys = {
            frame: role_for_confirmed_frame(frame, ramp_start, high_start, len(frames))
            for frame in confirmed_frames
        }
        key_source = "command_line" if cli_keys is not None else "registry"
    else:
        selected_keys = {}
        key_source = "not_confirmed"

    phase = np.full(len(frames), "high-temperature period", dtype=object)
    phase[:ramp_start] = "low-temperature period"
    phase[ramp_start:high_start] = "heating ramp"
    acquisition_time = temperature["acquisition_time_s"].to_numpy(float)
    ramp_time = acquisition_time[ramp_start]
    high_time = acquisition_time[high_start]

    manifest = pd.DataFrame(
        {
            "frame": frames,
            "raw_frame_number": srs["frameNo"].astype(int),
            "acquisition_time_s": acquisition_time,
            "measured_temperature_C": temperature["measured_temperature_C"],
            "nominal_stage_C": np.nan,
            "phase": phase,
            "time_since_heating_start_s": np.where(
                np.arange(len(frames)) >= ramp_start,
                acquisition_time - ramp_time,
                np.nan,
            ),
            "time_at_high_temperature_s": np.where(
                np.arange(len(frames)) >= high_start,
                acquisition_time - high_time,
                np.nan,
            ),
            "d5i_monitor": monitor,
            "d5i_normalization_factor": scale,
            "exposure_s": srs["count_time"].to_numpy(float),
            "transmission": srs["transmission"].to_numpy(float),
            "key_frame_suggested": [frame in suggested_keys for frame in frames],
            "key_frame_selected": [frame in selected_keys for frame in frames],
            "key_frame_role": [selected_keys.get(frame, "") for frame in frames],
            "key_frame_selection_source": key_source,
        }
    )

    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else PROJECT_ROOT
        / "data"
        / entry["sample_slug"]
        / f"scan_{args.scan}"
        / "processed"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    planned = [
        output_dir / f"scan{args.scan}_{cut}_{len(frames)}frames_{kind}.txt"
        for cut in CUTS
        for kind in ("raw", "norm")
    ]
    planned.extend(
        [
            output_dir / f"scan{args.scan}_frame_manifest.csv",
            output_dir / f"scan{args.scan}_input_provenance.json",
            output_dir / f"scan{args.scan}_key_frame_suggestions.csv",
        ]
    )
    existing = [path for path in planned if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing prepared files. Use --overwrite only "
            "after checking the target:\n  " + "\n  ".join(map(str, existing))
        )

    for cut in CUTS:
        q, intensity = profiles[cut]
        raw = pd.concat(
            [pd.Series(q, name="q_Ainv"), intensity.reset_index(drop=True)],
            axis=1,
        )
        normalized = raw.copy()
        for index, frame in enumerate(frames):
            normalized[str(frame)] = raw[str(frame)] * scale[index]
        raw_path = output_dir / f"scan{args.scan}_{cut}_{len(frames)}frames_raw.txt"
        norm_path = output_dir / f"scan{args.scan}_{cut}_{len(frames)}frames_norm.txt"
        raw.to_csv(raw_path, sep="\t", index=False)
        normalized.to_csv(norm_path, sep="\t", index=False)
        outputs[f"{cut}_raw"] = raw_path
        outputs[f"{cut}_norm"] = norm_path

    manifest_path = output_dir / f"scan{args.scan}_frame_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    outputs["manifest"] = manifest_path

    suggestion_rows = []
    for frame, role in suggested_keys.items():
        row = manifest.loc[manifest["frame"] == frame].iloc[0]
        suggestion_rows.append(
            {
                "frame": frame,
                "measured_temperature_C": row["measured_temperature_C"],
                "acquisition_time_s": row["acquisition_time_s"],
                "phase": row["phase"],
                "suggested_role": role,
                "selected_for_fitting": frame in selected_keys,
                "selection_source": key_source,
                "suggestion_source": suggestion_source,
            }
        )
    suggestions_path = output_dir / f"scan{args.scan}_key_frame_suggestions.csv"
    pd.DataFrame(suggestion_rows).to_csv(suggestions_path, index=False)
    outputs["key_frame_suggestions"] = suggestions_path

    q_step = float(np.median(np.diff(q_reference)))
    provenance = {
        "scan": int(args.scan),
        "sample_label": entry["sample_label"],
        "sample_slug": entry["sample_slug"],
        "frame_count": len(frames),
        "registry": str(registry_path),
        "sources": {
            "temperature_table": str(temperature_path),
            "raw_point_metadata": str(dat_path),
            "profiles": {
                cut: {
                    "path": str(profile_sources[cut]["path"]),
                    "sheet": profile_sources[cut]["sheet"],
                }
                for cut in CUTS
            },
        },
        "validation": {
            "zero_based_profile_frames": True,
            "raw_frame_mapping": "processed frame = SRS frameNo - 1",
            "maximum_temperature_alignment_difference_C": float(
                np.nanmax(temperature_difference)
            ),
            "identical_q_grid_across_cuts": True,
            "q_rows": int(len(q_reference)),
            "q_min_Ainv": float(q_reference[0]),
            "q_max_Ainv": float(q_reference[-1]),
            "q_step_Ainv": q_step,
            "finite_intensity_coverage_by_cut": profile_coverage,
            "missing_value_handling": "The fitter excludes non-finite q/intensity pairs within each fitting window.",
        },
        "intensity_normalization": {
            "formula": "I_norm(frame,q) = I_raw(frame,q) * median(d5i) / d5i(frame)",
            "d5i_median": monitor_median,
            "note": "Multiplicative frame normalization changes area scale but not q0 or FWHM within a frame.",
        },
        "phase_suggestion": {
            "ramp_start_frame": ramp_start,
            "high_temperature_start_frame": high_start,
            **phase_details,
            "warning": "Phase labels are temperature-trace suggestions and must be checked against the experimental protocol/log.",
        },
        "key_frames": {
            "automatic_temperature_suggestion": list(automatic_keys),
            "reviewed_registry_suggestion": list(suggested_keys),
            "suggestion_source": suggestion_source,
            "selected": list(selected_keys),
            "selection_source": key_source,
            "confirmed": key_source in {"registry", "command_line"},
        },
        "analysis_boundary": "Peak presence, q0/d, and integrated area. No CCL interpretation for secondary samples.",
    }
    provenance_path = output_dir / f"scan{args.scan}_input_provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    outputs["provenance"] = provenance_path
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate and prepare one secondary-sample GIWAXS scan."
    )
    parser.add_argument("--scan", required=True, help="Beamline scan number.")
    parser.add_argument(
        "--registry", type=Path, default=DEFAULT_REGISTRY, help="Scan source registry."
    )
    parser.add_argument(
        "--key-frames",
        type=parse_frame_list,
        default=None,
        help="Optional confirmed comma-separated key frames; overrides the registry.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional prepared-data destination.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing prepared outputs after explicit review.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    outputs = prepare(args)
    print("Prepared scan without fitting.")
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
