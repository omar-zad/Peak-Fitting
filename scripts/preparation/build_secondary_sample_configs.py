#!/usr/bin/env python3
"""Build candidate or profile-locked configs for a secondary GIWAXS scan.

Candidate configs inherit the deliberately limited peak-window library already
tested on scan 587214, but are labelled provisional and must be reviewed on
key-frame overlays before use. The ``--freeze-from-results`` mode copies the
profile family selected for each window in a completed key-frame run so the
same family and pseudo-Voigt eta are used for an optional all-frame run.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CUTS = ("FR", "IP", "OOP")
DEFAULT_TEMPLATE_DIR = (
    PROJECT_ROOT / "configs" / "in_situ" / "scan_587214" / "representative"
)


def parse_frames(value: str | None) -> list[int] | None:
    if value is None:
        return None
    frames = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not frames or len(frames) != len(set(frames)):
        raise argparse.ArgumentTypeError(
            "--frames requires a unique comma-separated frame list."
        )
    return frames


def load_manifest(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    required = {"frame", "measured_temperature_C", "phase"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Manifest is missing: {', '.join(missing)}")
    table["frame"] = pd.to_numeric(table["frame"], errors="raise").astype(int)
    table["measured_temperature_C"] = pd.to_numeric(
        table["measured_temperature_C"], errors="raise"
    ).astype(float)
    if table["frame"].tolist() != list(range(len(table))):
        raise ValueError("Manifest must contain continuous zero-based frames.")
    return table


def select_frames(manifest: pd.DataFrame, override: list[int] | None) -> list[int]:
    if override is not None:
        frames = override
    elif "key_frame_selected" in manifest:
        selected = manifest["key_frame_selected"].astype(str).str.lower().isin(
            {"true", "1", "yes"}
        )
        frames = manifest.loc[selected, "frame"].tolist()
    else:
        frames = []
    if not frames:
        raise ValueError(
            "No key frames are selected. Review the suggestion CSV and pass "
            "--frames, or update the registry and rebuild the prepared data."
        )
    invalid = sorted(set(frames) - set(manifest["frame"]))
    if invalid:
        raise ValueError(f"Selected frames are absent from the manifest: {invalid}")
    return frames


def high_anchor_frame(manifest: pd.DataFrame, frames: list[int]) -> int:
    selected = manifest.set_index("frame").loc[frames].copy()
    high = selected[
        selected["phase"].astype(str).str.contains("high", case=False, na=False)
    ]
    if high.empty:
        high = selected
    # Prefer a representative high-temperature frame rather than an extreme
    # final frame; ties are resolved toward the acquisition midpoint.
    maximum = float(high["measured_temperature_C"].max())
    near_max = high[high["measured_temperature_C"] >= maximum - 1.0]
    candidates = near_max.index.to_list() or high.index.to_list()
    return int(candidates[len(candidates) // 2])


def template_path(template_dir: Path, cut: str) -> Path:
    matches = sorted(template_dir.glob(f"*_{cut}_keyframes.json"))
    if len(matches) != 1:
        found = ", ".join(path.name for path in matches) or "none"
        raise FileNotFoundError(
            f"Expected one {cut} key-frame template in {template_dir}; found {found}."
        )
    return matches[0]


def candidate_config(
    source: Path,
    manifest_path: Path,
    manifest: pd.DataFrame,
    frames: list[int],
    scan: str,
    sample_label: str,
    cut: str,
) -> dict[str, Any]:
    config = json.loads(source.read_text(encoding="utf-8"))
    config["frame_temperatures_C"] = {
        str(int(row.frame)): float(row.measured_temperature_C)
        for row in manifest.itertuples(index=False)
    }
    config["representative_frames"] = frames
    config["analysis_scope"] = (
        "Candidate-window key-frame screening for peak presence, q0/d and "
        "integrated area. No CCL interpretation."
    )
    anchor = high_anchor_frame(manifest, frames)
    anchor_temperature = float(
        manifest.set_index("frame").loc[anchor, "measured_temperature_C"]
    )
    for window in config["windows"]:
        window.pop("profiles_to_compare", None)
        window.pop("eta_fixed", None)
        window["candidate_window_requires_visual_review"] = True
        for peak in window["peaks"]:
            peak["anchor_frame"] = anchor
            peak["anchor_temperature_C"] = anchor_temperature
            peak["provisional_assignment"] = True

    config["series_metadata"] = {
        "scan": int(scan),
        "sample_label": sample_label,
        "cut": cut,
        "series_type": "secondary_in_situ_time_series",
        "frame_manifest": str(manifest_path),
        "key_frames": frames,
        "anchor_frame": anchor,
        "anchor_temperature_C": anchor_temperature,
        "template_config": str(source),
        "scientific_status": (
            "Candidate windows transferred from scan 587214 for consistent "
            "screening; peak count, bounds and background must be accepted or "
            "revised from this sample's own key-frame overlays."
        ),
    }
    return config


def freeze_profiles(
    config: dict[str, Any],
    model_selection_path: Path,
) -> dict[str, Any]:
    selection = pd.read_csv(model_selection_path)
    required = {"window", "profile", "eta_shared", "selected"}
    missing = sorted(required - set(selection.columns))
    if missing:
        raise ValueError(
            f"{model_selection_path} is missing: {', '.join(missing)}"
        )
    selected_mask = selection["selected"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    selected = selection.loc[selected_mask].copy()
    by_window = {str(row.window): row for row in selected.itertuples(index=False)}
    expected = {window["key"] for window in config["windows"]}
    missing_windows = sorted(expected - set(by_window))
    if missing_windows:
        raise ValueError(
            "No selected profile was recorded for windows: "
            + ", ".join(missing_windows)
        )
    lock_summary = {}
    for window in config["windows"]:
        row = by_window[window["key"]]
        profile = str(row.profile)
        window["profiles_to_compare"] = [profile]
        window.pop("eta_fixed", None)
        eta = None
        if profile == "pseudo_voigt":
            eta = float(row.eta_shared)
            if not 0.0 <= eta <= 1.0:
                raise ValueError(f"Invalid eta for {window['key']}: {eta}")
            window["eta_fixed"] = eta
        lock_summary[window["key"]] = {"profile": profile, "eta": eta}
    config["profile_lock"] = {
        "source_model_selection": str(model_selection_path),
        "window_profiles": lock_summary,
        "note": (
            "Profile families were selected across the approved key frames "
            "and are fixed here for all-frame tracking."
        ),
    }
    return config


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build candidate or key-frame-profile-locked scan configs."
    )
    parser.add_argument("--scan", required=True)
    parser.add_argument("--sample-label", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--template-dir", type=Path, default=DEFAULT_TEMPLATE_DIR
    )
    parser.add_argument("--frames", type=parse_frames, default=None)
    parser.add_argument(
        "--freeze-from-results",
        type=Path,
        default=None,
        help="Completed key-frame results root containing FR/IP/OOP subfolders.",
    )
    args = parser.parse_args()

    manifest_path = args.manifest.expanduser().resolve()
    manifest = load_manifest(manifest_path)
    frames = select_frames(manifest, args.frames)
    template_dir = args.template_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for cut in CUTS:
        source = template_path(template_dir, cut)
        config = candidate_config(
            source,
            manifest_path,
            manifest,
            frames,
            args.scan,
            args.sample_label,
            cut,
        )
        suffix = "candidate_keyframes"
        if args.freeze_from_results is not None:
            model_selection = (
                args.freeze_from_results.expanduser().resolve()
                / cut
                / "model_selection.csv"
            )
            if not model_selection.is_file():
                raise FileNotFoundError(model_selection)
            config = freeze_profiles(config, model_selection)
            suffix = "locked_allframes"
        output = output_dir / f"scan{args.scan}_{cut}_{suffix}.json"
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
