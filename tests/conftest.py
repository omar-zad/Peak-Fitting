"""Shared helpers for the test suite."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
FITTER = SCRIPTS / "fitting" / "fit_giwaxs_series.py"
EXAMPLE_DATA = REPO / "example_data"
CONFIGS = REPO / "configs"


def run_script(script: Path, *arguments: str | Path) -> subprocess.CompletedProcess[str]:
    """Run a repository script with the test interpreter and fail loudly on error."""
    completed = subprocess.run(
        [sys.executable, str(script), *map(str, arguments)],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        pytest.fail(
            f"{script.relative_to(REPO)} failed with exit code {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
        )
    return completed


def fit(input_file: Path, config: Path, outdir: Path, *extra: str) -> Path:
    """Run the fitter with --bootstrap 0 and return the output directory."""
    run_script(FITTER, input_file, "--config", config, "--outdir", outdir, "--bootstrap", "0", *extra)
    return outdir
