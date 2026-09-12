#!/usr/bin/env python3
"""Build checkpoint-revised candidate-v3 configs for scan 587225.

Version 3 preserves every accepted candidate-v2 window and explicitly applies
the changes listed in ``FIRST_FRAME_DIAGNOSTIC_V2_QC.md``.  It does not fit
data.  Narrow features marked ``position_only`` may support peak presence and
q position, but their apparent FWHM is intentionally not a scientific result.
Broad shoulder components marked ``reference_like_component`` exist only to
prevent an asymmetric/broad base from biasing the adjacent sharp-peak centre.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import pandas as pd

import build_scan587225_candidate_v2 as v2


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = v2.DEFAULT_MANIFEST
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT
    / "configs"
    / "secondary_samples"
    / "scan_587225"
    / "candidate_v3"
)


def tagged_peak(*args: Any, **kwargs: Any) -> dict[str, Any]:
    position_only = bool(kwargs.pop("position_only", False))
    reference_like = bool(kwargs.pop("reference_like_component", False))
    result = v2.peak(*args, **kwargs)
    if position_only:
        result["position_only"] = True
        result["suppress_fwhm_interpretation"] = True
    if reference_like:
        result["reference_like_component"] = True
        result["exclude_from_peak_assignment"] = True
        result["suppress_fwhm_interpretation"] = True
    return result


def copied_window(windows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    return copy.deepcopy(next(item for item in windows if item["key"] == key))


def fr_windows() -> list[dict[str, Any]]:
    old = v2.fr_windows()
    early_low = v2.common_low_q("FR")
    early_low["key"] = "low_q_pair_early"
    early_low["title"] = "FR early low-q broad/0.331 pair"
    early_low["temperature_max_C"] = v2.EARLY_MAX_C
    early_low["stage2_revision"] = (
        "The pair passed frame 0 but the sharp component was absent and bound-limited at frame 40."
    )

    late_low = v2.window(
        "low_q_broad_late",
        "FR late single broad low-q component",
        0.225,
        0.405,
        [
            tagged_peak(
                "low_q_broad_late",
                "FR late broad low-q component near 0.318",
                0.318,
                0.295,
                0.338,
                0.080,
                0.040,
                0.140,
                40,
                84.1,
            )
        ],
        "Frame 40 is explained by one broad envelope; the v2 0.331 component was not detected and reached its width floor.",
        temperature_min_C=v2.LATE_MIN_C,
    )

    early_0436 = v2.window(
        "feature_0436_early",
        "FR early peak near 0.436",
        0.405,
        0.480,
        [
            tagged_peak(
                "feature_0436_early",
                "FR early peak near 0.436",
                0.436,
                0.423,
                0.450,
                0.018,
                0.009,
                0.050,
                0,
                31.3,
            )
        ],
        "The v2 frame-0 fit strongly supported 0.436 A^-1; the unsupported 0.486 component was removed from this window.",
        temperature_max_C=55.0,
    )

    transition_0486 = v2.window(
        "feature_0486_transition",
        "FR transition-only feature near 0.48",
        0.445,
        0.515,
        [
            tagged_peak(
                "feature_0486_transition",
                "FR transition-only feature near 0.48",
                0.480,
                0.465,
                0.501,
                0.030,
                0.012,
                0.070,
                27,
                63.9,
            )
        ],
        "The weak 0.486 candidate is tested only in selected ramp frames 24 and 27, not forced at frame 0.",
        temperature_min_C=45.0,
        temperature_max_C=v2.EARLY_MAX_C,
    )

    late_0455 = v2.window(
        "feature_0455_late",
        "FR late peak near 0.455",
        0.425,
        0.500,
        [
            tagged_peak(
                "feature_0455_late",
                "FR late peak near 0.455",
                0.455,
                0.443,
                0.467,
                0.018,
                0.009,
                0.050,
                40,
                84.1,
            )
        ],
        "Separated from 0.533 after the v2 joint window placed the non-detected 0.533 width at its ceiling.",
        temperature_min_C=v2.LATE_MIN_C,
    )

    late_0533 = v2.window(
        "feature_0533_late",
        "FR late broad candidate near 0.533",
        0.490,
        0.585,
        [
            tagged_peak(
                "feature_0533_late",
                "FR late broad candidate near 0.533",
                0.533,
                0.515,
                0.550,
                0.045,
                0.018,
                0.100,
                70,
                84.3,
            )
        ],
        "Independent window permits genuine non-detection at frame 40 without biasing 0.455; frame 70 is the configured reference.",
        temperature_min_C=v2.LATE_MIN_C,
    )

    revised_1538 = v2.window(
        "feature_1538_early_components",
        "FR early broad underlay and narrow peak near 1.538",
        1.470,
        1.590,
        [
            tagged_peak(
                "feature_1515_early_underlay",
                "FR broad underlay below the 1.538 apex",
                1.515,
                1.492,
                1.528,
                0.055,
                0.025,
                0.110,
                0,
                31.3,
                reference_like_component=True,
            ),
            tagged_peak(
                "feature_1538_early_sharp",
                "FR early narrow peak near 1.538 (position only)",
                1.538,
                1.529,
                1.549,
                0.015,
                0.006,
                0.035,
                0,
                31.3,
                position_only=True,
            )
        ],
        "The v3 quadratic-background checkpoint still followed the broad underlay and left a standardized residual of 10.57; an explicit reference-like underlay now protects the sharp-apex centre.",
        temperature_max_C=v2.EARLY_MAX_C,
    )

    revised_high_q = v2.window(
        "high_q_early_components",
        "FR early broad base, 1.999 apex and 2.08 peak",
        1.930,
        2.160,
        [
            tagged_peak(
                "feature_1980_early_base",
                "FR broad lower-q base below 1.999",
                1.978,
                1.955,
                1.991,
                0.055,
                0.025,
                0.110,
                0,
                31.3,
                reference_like_component=True,
            ),
            tagged_peak(
                "feature_1999_early_sharp",
                "FR sharp early apex near 1.999 (position only)",
                1.999,
                1.992,
                2.010,
                0.014,
                0.006,
                0.035,
                0,
                31.3,
                position_only=True,
            ),
            tagged_peak(
                "feature_2080_early",
                "FR early feature near 2.08",
                2.080,
                2.050,
                2.105,
                0.045,
                0.015,
                0.090,
                0,
                31.3,
            ),
        ],
        "The v2 two-component window left a standardized residual of 6.53 and biased the sharp-apex centroid to 1.990647 A^-1.",
        temperature_max_C=v2.EARLY_MAX_C,
    )

    return [
        early_low,
        late_low,
        early_0436,
        transition_0486,
        late_0455,
        late_0533,
        copied_window(old, "feature_0754_early"),
        copied_window(old, "feature_0941_early_broad"),
        revised_1538,
        copied_window(old, "feature_1720_early"),
        revised_high_q,
    ]


def ip_windows() -> list[dict[str, Any]]:
    old = v2.ip_windows()
    revised_1544 = v2.window(
        "feature_1544_early_sharp",
        "IP early narrow peak near 1.544 over curved local underlay",
        1.470,
        1.590,
        [
            tagged_peak(
                "feature_1544_early_sharp",
                "IP early narrow peak near 1.544 (position only)",
                1.544,
                1.531,
                1.558,
                0.015,
                0.006,
                0.035,
                0,
                31.3,
                position_only=True,
            )
        ],
        "A quadratic local background addresses the v2 structured residual from the broad lower-q underlay.",
        temperature_max_C=v2.EARLY_MAX_C,
    )
    revised_1544["background_order"] = 2
    return [
        copied_window(old, "low_q_pair"),
        copied_window(old, "mid_q_moving"),
        revised_1544,
        copied_window(old, "feature_1720_early"),
        copied_window(old, "feature_2080_early"),
    ]


def oop_windows() -> list[dict[str, Any]]:
    old = v2.oop_windows()
    low_broad = v2.window(
        "low_q_broad",
        "OOP single broad low-q component",
        0.225,
        0.405,
        [
            tagged_peak(
                "low_q_broad",
                "OOP broad low-q component near 0.32",
                0.320,
                0.295,
                0.338,
                0.085,
                0.040,
                0.140,
                40,
                84.1,
            )
        ],
        "Both v2 low-q components were centre-bound and not detected at frames 0 and 40; one broad envelope is supported.",
    )

    feature_0835 = copied_window(old, "feature_0835")
    narrow_0835 = feature_0835["peaks"][0]
    narrow_0835["fwhm_guess"] = 0.012
    narrow_0835["fwhm_min"] = 0.006
    narrow_0835["fwhm_max"] = 0.050
    narrow_0835["position_only"] = True
    narrow_0835["suppress_fwhm_interpretation"] = True
    feature_0835["stage2_revision"] = (
        "One-q-step lower width bound permits centre estimation; FWHM remains under-sampled and unreported."
    )

    early_high = copied_window(old, "high_q_early_pair")
    for component in early_high["peaks"]:
        component["fwhm_min"] = 0.006
        component["fwhm_guess"] = min(component["fwhm_guess"], 0.018)
        component["fwhm_max"] = min(component["fwhm_max"], 0.060)
        component["position_only"] = True
        component["suppress_fwhm_interpretation"] = True
    early_high["stage2_revision"] = (
        "The 1.75 component was visibly present but width-bound at 0.009 A^-1; one-q-step width is allowed only for q-position/presence screening."
    )

    revised_2223 = v2.window(
        "feature_2223_early_components",
        "OOP early broad shoulder and sharp 2.223 apex",
        2.145,
        2.290,
        [
            tagged_peak(
                "feature_2190_early_shoulder",
                "OOP broad lower-q shoulder below 2.223",
                2.190,
                2.168,
                2.206,
                0.045,
                0.020,
                0.100,
                0,
                31.3,
                reference_like_component=True,
            ),
            tagged_peak(
                "feature_2223_early_sharp",
                "OOP sharp early apex near 2.223 (position only)",
                2.223,
                2.207,
                2.240,
                0.018,
                0.006,
                0.050,
                0,
                31.3,
                position_only=True,
            ),
        ],
        "The v2 single symmetric line left reduced chi-square 11.86 and a standardized residual of 7.02; the raw profile contains a lower-q shoulder.",
        temperature_max_C=v2.EARLY_MAX_C,
    )

    return [
        low_broad,
        copied_window(old, "feature_0436_early"),
        copied_window(old, "feature_0536_late"),
        copied_window(old, "feature_0754_early"),
        feature_0835,
        copied_window(old, "feature_1495_early"),
        early_high,
        revised_2223,
    ]


def update_metadata(config: dict[str, Any], cut: str) -> None:
    metadata = config["series_metadata"]
    metadata["candidate_version"] = 3
    metadata["source_of_windows"] = (
        "Scan-587225 Stage-1 visual review plus candidate-v2 frame-0/frame-40 numerical and overlay QC."
    )
    metadata["parent_candidate"] = str(
        PROJECT_ROOT
        / "configs"
        / "secondary_samples"
        / "scan_587225"
        / "candidate_v2"
        / f"scan587225_{cut}_candidate_keyframes.json"
    )
    metadata["parent_qc_report"] = str(
        PROJECT_ROOT
        / "results"
        / "dimitar_n2n2"
        / "scan_587225"
        / "FIRST_FRAME_DIAGNOSTIC_V2_QC.md"
    )
    metadata["width_policy"] = (
        "Features marked position_only or suppress_fwhm_interpretation may support presence/q only. Apparent FWHM is diagnostic and must not be used as a result."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build checkpoint-revised candidate-v3 configs for scan 587225."
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    manifest_path = args.manifest.expanduser().resolve()
    manifest = pd.read_csv(manifest_path)
    selected = manifest.loc[
        manifest["key_frame_selected"].astype(str).str.lower().isin(
            {"true", "1", "yes"}
        ),
        "frame",
    ].astype(int).tolist()
    if selected != v2.SELECTED_FRAMES:
        raise ValueError(
            f"Manifest key frames {selected} do not match {v2.SELECTED_FRAMES}."
        )

    window_sets = {"FR": fr_windows(), "IP": ip_windows(), "OOP": oop_windows()}
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for cut, windows in window_sets.items():
        config = v2.build_config(cut, windows, manifest, manifest_path)
        config["analysis_scope"] = (
            "Checkpoint-revised scan-specific candidate-v3 key-frame screening "
            "for peak presence and q0/d. Under-resolved/reference-like widths "
            "are suppressed; no CCL interpretation."
        )
        update_metadata(config, cut)
        output = output_dir / f"scan587225_{cut}_candidate_keyframes.json"
        output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
