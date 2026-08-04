#!/usr/bin/env python3
"""Apply an auditable manual visual-QC exclusion list to accepted q positions."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KEY_COLUMNS = ["cut", "frame", "window", "peak"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    args = parser.parse_args()

    report_dir = args.report_dir.expanduser().resolve()
    accepted_path = report_dir / "accepted_peak_positions.csv"
    exclusions_path = args.exclusions.expanduser().resolve()
    accepted = pd.read_csv(accepted_path)
    exclusions = pd.read_csv(exclusions_path)
    missing = sorted(set(KEY_COLUMNS + ["reason"]) - set(exclusions.columns))
    if missing:
        raise ValueError(f"{exclusions_path} is missing: {', '.join(missing)}")

    exclusions = exclusions.loc[:, KEY_COLUMNS + ["reason"]].copy()
    exclusions["frame"] = pd.to_numeric(exclusions["frame"], errors="raise").astype(int)
    accepted["frame"] = pd.to_numeric(accepted["frame"], errors="raise").astype(int)
    marked = accepted.merge(
        exclusions,
        on=KEY_COLUMNS,
        how="left",
        validate="one_to_one",
    )
    unmatched = exclusions.merge(
        accepted.loc[:, KEY_COLUMNS],
        on=KEY_COLUMNS,
        how="left",
        indicator=True,
    ).loc[lambda table: table["_merge"] == "left_only"]
    if not unmatched.empty:
        raise ValueError(
            "Manual exclusion did not match an accepted position row:\n"
            + unmatched.loc[:, KEY_COLUMNS].to_string(index=False)
        )

    excluded = marked.loc[marked["reason"].notna()].copy()
    retained = marked.loc[marked["reason"].isna()].drop(columns=["reason"])
    retained.to_csv(report_dir / "accepted_peak_positions_after_visual_qc.csv", index=False)
    excluded.to_csv(report_dir / "manual_qc_excluded_positions.csv", index=False)
    print(
        f"Retained {len(retained)} accepted positions; manually excluded "
        f"{len(excluded)} in {report_dir}"
    )


if __name__ == "__main__":
    main()
