"""Regression test: the committed Drop40 FR checkpoint must be reproducible.

The repository ships a verified reference run under
``example_results/drop40_FR_checkpoint/``. This test re-runs the exact command
from the README against the committed example data and asserts that the
scientific decisions and the fitted numbers still match.

Tolerances
----------
The reference was produced on macOS (arm64, Python 3.13, numpy 2.5.1). Re-runs
on Linux under Python 3.10 (numpy 2.2.6 / pandas 2.3.3 / scipy 1.15.3) and
Python 3.13 (numpy 2.4.6 / pandas 3.0.5 / scipy 1.17.1) reproduced it with a
maximum absolute difference of 1.1e-9 in fitted peak centres and 6.7e-11 in
BIC. The tolerances below sit well above that observed spread but far below any
difference that could change a scientific conclusion.

Near-tie windows
----------------
Three windows in this checkpoint are decided by a BIC margin below one unit:
``pi_pi_sharp_pair`` (-0.55, where pseudo-Voigt has the lower raw BIC but the
documented parsimony rule retains the simpler Gaussian), ``q0754`` (+0.70) and
``q1999_single`` (+0.58). Observed cross-platform BIC variation is about nine
orders of magnitude smaller, so exact assertions remain safe. A failure on one
of those windows indicates a real near-tie flip worth investigating, not
numerical noise.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
REFERENCE_DIR = REPO / "example_results" / "drop40_FR_checkpoint"
INPUT_FILE = REPO / "example_data" / "drop40" / "drop40_FR_9frames_norm.txt"
CONFIG_FILE = REPO / "configs" / "drop40" / "peakfit_config_FR_revised.json"
FITTER = REPO / "scripts" / "fitting" / "fit_giwaxs_series.py"

Q_ABSOLUTE_TOLERANCE = 1e-6
FWHM_ABSOLUTE_TOLERANCE = 1e-5
AREA_RELATIVE_TOLERANCE = 1e-4
BIC_ABSOLUTE_TOLERANCE = 1e-6

MERGE_KEYS = ["window", "peak"]


@pytest.fixture(scope="module")
def rerun_directory(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Run the README's first-frame checkpoint into a temporary directory."""
    output_dir = tmp_path_factory.mktemp("drop40_FR_checkpoint")
    completed = subprocess.run(
        [
            sys.executable,
            str(FITTER),
            str(INPUT_FILE),
            "--config",
            str(CONFIG_FILE),
            "--outdir",
            str(output_dir),
            "--stage",
            "first",
            "--bootstrap",
            "0",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        pytest.fail(
            "The checkpoint run failed.\n"
            f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
        )
    return output_dir


def _merged(rerun_directory: Path, filename: str) -> pd.DataFrame:
    reference = pd.read_csv(REFERENCE_DIR / filename)
    rerun = pd.read_csv(rerun_directory / filename)
    keys = MERGE_KEYS if filename == "peak_parameters.csv" else ["window", "profile"]
    merged = reference.merge(rerun, on=keys, suffixes=("_reference", "_rerun"))
    assert len(merged) == len(reference), (
        f"{filename}: the re-run produced a different set of rows "
        f"({len(merged)} matched out of {len(reference)} reference rows)"
    )
    return merged


def test_example_files_are_present() -> None:
    for path in (REFERENCE_DIR, INPUT_FILE, CONFIG_FILE, FITTER):
        assert path.exists(), f"missing repository file: {path}"


def test_selected_profile_family_matches_reference(rerun_directory: Path) -> None:
    """The per-window line-shape decision is the primary scientific output."""
    merged = _merged(rerun_directory, "model_selection.csv")
    selected = merged[merged.selected_reference | merged.selected_rerun]
    disagreements = selected[selected.selected_reference != selected.selected_rerun]
    assert disagreements.empty, (
        "profile-family selection changed for: "
        + ", ".join(sorted(disagreements.window.unique()))
    )


def test_bic_values_match_reference(rerun_directory: Path) -> None:
    merged = _merged(rerun_directory, "model_selection.csv")
    difference = (merged.bic_rerun - merged.bic_reference).abs()
    worst = difference.max()
    assert worst < BIC_ABSOLUTE_TOLERANCE, (
        f"largest BIC difference {worst:.3e} exceeds {BIC_ABSOLUTE_TOLERANCE:.0e} "
        f"in window {merged.loc[difference.idxmax(), 'window']}"
    )


def test_peak_centres_match_reference(rerun_directory: Path) -> None:
    merged = _merged(rerun_directory, "peak_parameters.csv")
    difference = (merged.q0_fit_Ainv_rerun - merged.q0_fit_Ainv_reference).abs()
    worst = difference.max()
    assert worst < Q_ABSOLUTE_TOLERANCE, (
        f"largest fitted-centre difference {worst:.3e} 1/A exceeds "
        f"{Q_ABSOLUTE_TOLERANCE:.0e} for peak "
        f"{merged.loc[difference.idxmax(), 'peak']}"
    )


def test_peak_widths_and_areas_match_reference(rerun_directory: Path) -> None:
    merged = _merged(rerun_directory, "peak_parameters.csv")
    width_difference = (
        merged.apparent_fwhm_fit_Ainv_rerun - merged.apparent_fwhm_fit_Ainv_reference
    ).abs()
    assert width_difference.max() < FWHM_ABSOLUTE_TOLERANCE, (
        f"largest apparent-FWHM difference {width_difference.max():.3e} 1/A "
        f"exceeds {FWHM_ABSOLUTE_TOLERANCE:.0e}"
    )
    relative_area_difference = (
        (merged.area_total_rerun - merged.area_total_reference).abs()
        / merged.area_total_reference.abs()
    )
    assert relative_area_difference.max() < AREA_RELATIVE_TOLERANCE, (
        f"largest relative area difference {relative_area_difference.max():.3e} "
        f"exceeds {AREA_RELATIVE_TOLERANCE:.0e}"
    )


def test_detection_status_matches_reference(rerun_directory: Path) -> None:
    """Detection status drives what may be claimed, so it must be stable."""
    merged = _merged(rerun_directory, "peak_parameters.csv")
    disagreements = merged[
        merged.detection_status_reference != merged.detection_status_rerun
    ]
    assert disagreements.empty, (
        "detection status changed for: "
        + ", ".join(
            f"{row.window}/{row.peak} "
            f"({row.detection_status_reference} -> {row.detection_status_rerun})"
            for row in disagreements.itertuples()
        )
    )


def test_reported_positions_pass_the_same_gate(rerun_directory: Path) -> None:
    """The conservative reporting gate must admit exactly the same peaks."""
    merged = _merged(rerun_directory, "peak_parameters.csv")
    reference_reported = set(
        merged.loc[merged.q0_reported_Ainv_reference.notna()]
        .set_index(MERGE_KEYS)
        .index
    )
    rerun_reported = set(
        merged.loc[merged.q0_reported_Ainv_rerun.notna()].set_index(MERGE_KEYS).index
    )
    assert reference_reported == rerun_reported, (
        "the reporting gate admitted a different set of peaks; "
        f"only in reference: {sorted(reference_reported - rerun_reported)}; "
        f"only in re-run: {sorted(rerun_reported - reference_reported)}"
    )
