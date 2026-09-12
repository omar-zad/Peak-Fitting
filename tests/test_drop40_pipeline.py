"""Run the documented Drop40 chain end to end and check the numbers quoted in the README.

This module is marked slow: it fits all nine frames for FR, IP and OOP, the IP
one- and three-component alternatives and the two cross-sector sensitivity
configurations, then builds the QC package and both sensitivity summaries.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from conftest import CONFIGS, EXAMPLE_DATA, REPO, SCRIPTS, fit, run_script

pytestmark = pytest.mark.slow

DROP40 = EXAMPLE_DATA / "drop40"
CROSS_SECTOR_FRAMES = "197,198,199,200"
README_BIC = {"one": 905.2, "two": 474.8, "three": 503.2}


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("drop40")
    for cut in ("FR", "IP", "OOP"):
        fit(DROP40 / f"drop40_{cut}_9frames_norm.txt", CONFIGS / "drop40" / f"peakfit_config_{cut}_revised.json",
            root / f"drop40_{cut}_all", "--stage", "all", "--confirm-first-fit")
    for n in ("1", "3"):
        fit(DROP40 / "drop40_IP_9frames_norm.txt", CONFIGS / "drop40" / "sensitivity" / f"peakfit_config_IP_{n}component.json",
            root / f"IP_{n}component", "--stage", "all", "--confirm-first-fit")
    for cut, name in (("IP", "peakfit_config_IP_q0480_cross_sector.json"), ("OOP", "peakfit_config_OOP_q0459_cross_sector.json")):
        fit(DROP40 / f"drop40_{cut}_9frames_norm.txt", CONFIGS / "drop40" / "sensitivity" / name,
            root / "cross_sector_0459_0480_sensitivity" / cut, "--stage", "selected", "--frames", CROSS_SECTOR_FRAMES)
    return root


def test_ip_component_count_matches_readme_and_committed_comparison(runs: Path) -> None:
    outdir = runs / "IP_component_count_comparison"
    run_script(
        SCRIPTS / "reporting" / "compare_component_models.py",
        f"--run=one={runs / 'IP_1component'}", f"--run=two={runs / 'drop40_IP_all'}", f"--run=three={runs / 'IP_3component'}",
        "--outdir", outdir, "--title", "Drop40 IP low-q component-count decision",
    )
    table = pd.read_csv(outdir / "component_model_comparison.csv").set_index("label")
    for label, quoted in README_BIC.items():
        assert abs(table.loc[label, "bic"] - quoted) < 0.05, f"{label}-component BIC {table.loc[label, 'bic']:.2f} != README {quoted}"
    assert table["preferred"].idxmax() == "two"
    committed = pd.read_csv(REPO / "example_results" / "drop40_IP_component_count" / "component_model_comparison.csv").set_index("label")
    assert (table["bic"] - committed["bic"]).abs().max() < 1e-6


def test_qc_package_builds_from_the_nine_frame_runs(runs: Path) -> None:
    outdir = runs / "qc_report"
    run_script(
        SCRIPTS / "reporting" / "build_nine_frame_qc_report.py",
        *[f"--input={cut}={runs / f'drop40_{cut}_all'}" for cut in ("FR", "IP", "OOP")],
        "--output", outdir,
    )
    assert (outdir / "QC_REPORT.md").is_file()
    gates = pd.read_csv(outdir / "reporting_gate_summary.csv")
    assert not gates.empty
    combined = pd.read_csv(outdir / "peak_parameters_combined_with_gates.csv")
    assert set(combined["cut"]) == {"FR", "IP", "OOP"}


def test_cross_sector_summary_builds(runs: Path) -> None:
    root = runs / "cross_sector_0459_0480_sensitivity"
    run_script(SCRIPTS / "reporting" / "summarize_drop40_cross_sector_sensitivity.py", "--root", root)
    summary = pd.read_csv(root / "cross_sector_sensitivity_summary.csv")
    assert len(summary) == 8 and set(summary["cut"]) == {"IP", "OOP"}
    assert (root / "CROSS_SECTOR_SENSITIVITY_REPORT.md").is_file()
