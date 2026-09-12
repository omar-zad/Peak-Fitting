"""The configuration builders must run from the repository and reproduce what is committed."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from conftest import CONFIGS, EXAMPLE_DATA, SCRIPTS, run_script

DROP40_BUILDER = SCRIPTS / "preparation" / "build_peakfit_cut_configs.py"
IN_SITU_BUILDER = SCRIPTS / "preparation" / "build_in_situ_scan_configs.py"
SECONDARY_BUILDER = SCRIPTS / "preparation" / "build_secondary_sample_configs.py"
EXAMPLE_MANIFEST = EXAMPLE_DATA / "secondary_scan" / "scan587178_frame_manifest.csv"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_drop40_builder_regenerates_every_committed_configuration(tmp_path: Path) -> None:
    """Running the builder into a scratch folder must reproduce configs/drop40 exactly."""
    run_script(DROP40_BUILDER, "--output-root", tmp_path)
    generated = sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*.json"))
    assert len(generated) == 15, generated
    for relative in generated:
        committed = CONFIGS / "drop40" / relative
        assert committed.is_file(), f"builder wrote {relative} but it is not committed"
        assert _load(tmp_path / relative) == _load(committed), f"{relative} differs from the committed file"


def test_in_situ_builder_adapts_temperature_map_and_anchors(tmp_path: Path) -> None:
    run_script(IN_SITU_BUILDER, "--scan", "587178", "--manifest", EXAMPLE_MANIFEST, "--output-dir", tmp_path)
    manifest = pd.read_csv(EXAMPLE_MANIFEST)
    expected_map = {str(int(f)): float(t) for f, t in zip(manifest["frame"], manifest["measured_temperature_C"])}
    low, high = manifest["measured_temperature_C"].min(), manifest["measured_temperature_C"].max()
    for cut in ("FR", "IP", "OOP"):
        config = _load(tmp_path / f"scan587178_{cut}.json")
        assert config["frame_temperatures_C"] == expected_map
        anchors = {peak["anchor_temperature_C"] for window in config["windows"] for peak in window["peaks"]}
        assert anchors <= {low, high}, f"{cut}: anchors were not remapped to the measured extremes: {anchors}"


def test_secondary_builder_writes_reviewable_candidates(tmp_path: Path) -> None:
    frames = [0, 5, 10, 15, 20, 25, 60, 110]
    run_script(
        SECONDARY_BUILDER,
        "--scan", "587178", "--sample-label", "Example bc", "--manifest", EXAMPLE_MANIFEST,
        "--frames", ",".join(map(str, frames)), "--output-dir", tmp_path,
    )
    for cut in ("FR", "IP", "OOP"):
        config = _load(tmp_path / f"scan587178_{cut}_candidate_keyframes.json")
        assert config["representative_frames"] == frames
        assert config["series_metadata"]["key_frames"] == frames
        assert all(window["candidate_window_requires_visual_review"] for window in config["windows"])
        assert all(peak["provisional_assignment"] for window in config["windows"] for peak in window["peaks"])
        template = Path(config["series_metadata"]["template_config"])
        assert template.parent == CONFIGS / "secondary_samples" / "scan_587214" / "representative"
