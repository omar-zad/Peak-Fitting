"""Fast smoke tests for the reporting tools on first-frame fits of the example data."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from conftest import CONFIGS, EXAMPLE_DATA, SCRIPTS, fit, run_script

IP_INPUT = EXAMPLE_DATA / "drop40" / "drop40_IP_9frames_norm.txt"
IP_CONFIGS = {
    "one": CONFIGS / "drop40" / "sensitivity" / "peakfit_config_IP_1component.json",
    "two": CONFIGS / "drop40" / "peakfit_config_IP_revised.json",
    "three": CONFIGS / "drop40" / "sensitivity" / "peakfit_config_IP_3component.json",
}


@pytest.fixture(scope="module")
def first_frame_ip_runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("ip_component_runs")
    return {label: fit(IP_INPUT, config, root / label, "--stage", "first") for label, config in IP_CONFIGS.items()}


def test_component_comparison_checks_matching_windows(first_frame_ip_runs: dict[str, Path], tmp_path: Path) -> None:
    run_script(
        SCRIPTS / "reporting" / "compare_component_models.py",
        *[f"--run={label}={path}" for label, path in first_frame_ip_runs.items()],
        "--outdir", tmp_path, "--title", "smoke test",
    )
    table = pd.read_csv(tmp_path / "component_model_comparison.csv")
    assert list(table["label"]) == ["one", "two", "three"]
    assert list(table["n_components"]) == [1, 2, 3]
    assert table["n_points"].nunique() == 1, "all candidates must be compared on the same q points"
    assert table["preferred"].sum() == 1
    assert (tmp_path / "component_model_comparison.png").is_file()
