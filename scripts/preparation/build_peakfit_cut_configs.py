#!/usr/bin/env python3
"""Create explicit, versioned Drop40 fitting configurations for each cut."""

from __future__ import annotations

import copy
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_ROOT = PROJECT_ROOT / "configs"
BASE_CONFIG = CONFIG_ROOT / "base" / "peakfit_config.json"


def write_config(name: str, config: dict) -> None:
    if name.endswith("_revised.json"):
        directory = CONFIG_ROOT / "production"
    elif name == "peakfit_config_OOP.json":
        directory = CONFIG_ROOT / "base"
    else:
        directory = CONFIG_ROOT / "sensitivity"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {path}")


def lamellar_window(config: dict) -> dict:
    return next(
        window
        for window in config["windows"]
        if window["key"] == "lamellar_low_q"
    )


def window_by_key(config: dict, key: str) -> dict:
    return next(window for window in config["windows"] if window["key"] == key)


def high_temperature_pi_pi_window(two_components: bool = False) -> dict:
    """Describe the broad high-temperature pi-pi-region envelope."""
    peaks = [
        {
            "key": "q1710_broad_highT",
            "label": "high-temperature broad pi-pi-region envelope",
            "q_guess": 1.71,
            "q_min": 1.62,
            "q_max": 1.80,
            "fwhm_guess": 0.14,
            "fwhm_min": 0.05,
            "fwhm_max": 0.25,
            "anchor_temperature_C": 150,
            "provisional_assignment": True,
        }
    ]
    key = "pi_pi_highT_one_broad"
    title = "High-temperature broad pi-pi-region envelope"
    if two_components:
        key = "pi_pi_highT_two_broad"
        title = "High-temperature two-component pi-pi-region sensitivity"
        peaks = [
            {
                "key": "q1660_broad_highT",
                "label": "high-temperature lower-q broad component",
                "q_guess": 1.66,
                "q_min": 1.62,
                "q_max": 1.69,
                "fwhm_guess": 0.12,
                "fwhm_min": 0.05,
                "fwhm_max": 0.25,
                "anchor_temperature_C": 150,
                "provisional_assignment": True,
            },
            {
                "key": "q1750_broad_highT",
                "label": "high-temperature higher-q broad component",
                "q_guess": 1.75,
                "q_min": 1.70,
                "q_max": 1.80,
                "fwhm_guess": 0.12,
                "fwhm_min": 0.05,
                "fwhm_max": 0.25,
                "anchor_temperature_C": 150,
                "provisional_assignment": True,
            },
        ]
    return {
        "key": key,
        "title": title,
        "q_min": 1.60,
        "q_max": 1.90,
        "temperature_min_C": 65,
        "peaks": peaks,
    }


def revised_windows(config: dict) -> list[dict]:
    """Return the less-degenerate, scientifically scoped window set.

    The original broad windows combined peaks that appear at different
    temperatures and/or are separated by enough baseline to be fitted
    independently.  Splitting them reduces component trading without adding
    any temperature-specific constraints to the fitter.
    """
    heated_peaks = window_by_key(config, "heated_low_q")["peaks"]
    pi_pi_peaks = window_by_key(config, "pi_pi_cluster")["peaks"]
    high_q_peaks = window_by_key(config, "high_q_cluster")["peaks"]
    q2279_peak = copy.deepcopy(window_by_key(config, "q2279")["peaks"][0])

    replacements = {
        "heated_low_q": [
            {
                "key": "heated_low_q_0480_0536",
                "title": "Heating-induced q = 0.480 / 0.536 A^-1 pair",
                "q_min": 0.44,
                "q_max": 0.57,
                "peaks": copy.deepcopy(heated_peaks[:2]),
            },
            {
                "key": "heated_low_q_0604_0642",
                "title": "Heating-induced q = 0.604 / 0.642 A^-1 pair",
                "q_min": 0.57,
                "q_max": 0.69,
                "peaks": copy.deepcopy(heated_peaks[2:]),
            },
        ],
        "pi_pi_cluster": [
            {
                "key": "pi_pi_sharp_pair",
                "title": "Low-temperature sharp q = 1.713 / 1.750 A^-1 pair",
                "q_min": 1.66,
                "q_max": 1.82,
                "temperature_max_C": 60,
                "peaks": copy.deepcopy(pi_pi_peaks[1:3]),
            },
            high_temperature_pi_pi_window(two_components=True),
        ],
        "high_q_cluster": [
            {
                "key": "q1999_single",
                "title": "Nominal q = 1.999 A^-1 reflection",
                "q_min": 1.94,
                "q_max": 2.04,
                "peaks": [copy.deepcopy(high_q_peaks[0])],
            },
            {
                "key": "q2086_q2117_pair",
                "title": "Nominal q = 2.086 / 2.117 A^-1 pair",
                "q_min": 2.04,
                "q_max": 2.17,
                "peaks": copy.deepcopy(high_q_peaks[1:]),
            },
        ],
        "q2279": [
            {
                "key": "q2279",
                "title": "Early-temperature q = 2.279 A^-1 candidate",
                "q_min": 2.22,
                "q_max": 2.32,
                "peaks": [q2279_peak],
            }
        ],
    }

    revised: list[dict] = []
    for window in config["windows"]:
        replacement = replacements.get(window["key"])
        if replacement is None:
            revised.append(copy.deepcopy(window))
        else:
            revised.extend(replacement)
    return revised


def revised_ip_windows(config: dict) -> list[dict]:
    """Return an IP-specific window set without forcing FR/OOP components.

    IP retains the supported two-component lamellar model, but the remaining
    regions are narrowed or split so that peaks absent from a frame can be
    classified as such instead of trading area/width with unrelated broad
    components.
    """
    lamellar = copy.deepcopy(window_by_key(config, "lamellar_low_q"))
    lamellar["title"] = "IP low-q two-component candidate"
    lamellar["profile_note"] = (
        "Gaussian, Lorentzian, and pseudo-Voigt remain candidates; select one "
        "family across all nine IP frames by BIC."
    )

    q0754 = copy.deepcopy(window_by_key(config, "q0754"))
    q0754["title"] = "Early-temperature q = 0.754 A^-1 candidate"
    q0754["peaks"][0]["provisional_assignment"] = True

    q0947 = copy.deepcopy(window_by_key(config, "q0947"))
    q0947["title"] = "Conditional q = 0.947 A^-1 candidate"
    q0947["peaks"][0]["provisional_assignment"] = True

    q1532 = copy.deepcopy(window_by_key(config, "q1532"))
    q1532["title"] = "IP q approximately 1.54 A^-1 (position/area only)"
    q1532["reporting_scope"] = "q_position_and_area_only"
    q1532_peak = q1532["peaks"][0]
    q1532_peak.update(
        {
            "label": "IP approximately 1.54 (q/area only)",
            "q_guess": 1.540,
            "q_min": 1.515,
            "q_max": 1.555,
            "provisional_assignment": True,
        }
    )

    return [
        lamellar,
        {
            "key": "heated_low_q_ip",
            "title": "IP heating-induced low-q candidates",
            "q_min": 0.44,
            "q_max": 0.57,
            "peaks": [
                {
                    "key": "q0461_ip",
                    "label": "IP provisional approximately 0.461",
                    "q_guess": 0.461,
                    "q_min": 0.450,
                    "q_max": 0.480,
                    "fwhm_guess": 0.030,
                    "fwhm_min": 0.010,
                    "fwhm_max": 0.070,
                    "anchor_temperature_C": 150,
                    "provisional_assignment": True,
                },
                {
                    "key": "q0520_ip",
                    "label": "IP provisional approximately 0.520",
                    "q_guess": 0.520,
                    "q_min": 0.500,
                    "q_max": 0.545,
                    "fwhm_guess": 0.030,
                    "fwhm_min": 0.010,
                    "fwhm_max": 0.070,
                    "anchor_temperature_C": 150,
                    "provisional_assignment": True,
                },
            ],
        },
        q0754,
        q0947,
        q1532,
        {
            "key": "pi_pi_sharp_pair_ip",
            "title": "IP low-temperature sharp pi-pi-region pair",
            "q_min": 1.67,
            "q_max": 1.80,
            "temperature_max_C": 60,
            "peaks": [
                {
                    "key": "q1713_ip",
                    "label": "IP sharp approximately 1.713",
                    "q_guess": 1.713,
                    "q_min": 1.695,
                    "q_max": 1.730,
                    "fwhm_guess": 0.020,
                    "fwhm_min": 0.009,
                    "fwhm_max": 0.050,
                    "anchor_temperature_C": 40,
                    "provisional_assignment": True,
                },
                {
                    "key": "q1750_ip",
                    "label": "IP sharp approximately 1.750",
                    "q_guess": 1.750,
                    "q_min": 1.735,
                    "q_max": 1.765,
                    "fwhm_guess": 0.022,
                    "fwhm_min": 0.009,
                    "fwhm_max": 0.050,
                    "anchor_temperature_C": 40,
                    "provisional_assignment": True,
                },
            ],
        },
        high_temperature_pi_pi_window(two_components=True),
        {
            "key": "q2015_ip_provisional",
            "title": "IP provisional/unassigned q approximately 2.015 A^-1",
            "q_min": 1.96,
            "q_max": 2.05,
            "peaks": [
                {
                    "key": "q2015_ip",
                    "label": "IP provisional/unassigned approximately 2.015",
                    "q_guess": 2.015,
                    "q_min": 1.990,
                    "q_max": 2.035,
                    "fwhm_guess": 0.025,
                    "fwhm_min": 0.009,
                    "fwhm_max": 0.065,
                    "anchor_temperature_C": 40,
                    "provisional_assignment": True,
                }
            ],
        },
        {
            "key": "q2074_q2111_ip",
            "title": "IP provisional high-q sharp pair",
            "q_min": 2.04,
            "q_max": 2.15,
            "peaks": [
                {
                    "key": "q2074_ip",
                    "label": "IP provisional approximately 2.074",
                    "q_guess": 2.074,
                    "q_min": 2.055,
                    "q_max": 2.098,
                    "fwhm_guess": 0.025,
                    "fwhm_min": 0.009,
                    "fwhm_max": 0.070,
                    "anchor_temperature_C": 40,
                    "provisional_assignment": True,
                },
                {
                    "key": "q2111_ip",
                    "label": "IP provisional approximately 2.111",
                    "q_guess": 2.111,
                    "q_min": 2.100,
                    "q_max": 2.145,
                    "fwhm_guess": 0.022,
                    "fwhm_min": 0.009,
                    "fwhm_max": 0.060,
                    "anchor_temperature_C": 40,
                    "provisional_assignment": True,
                },
            ],
        },
    ]


def main() -> None:
    base = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))

    fr_two = copy.deepcopy(base)
    fr_two["analysis_variant"] = (
        "FR two-component low-q reference; norm-corrected nine-frame series"
    )
    write_config("peakfit_config_FR_2component.json", fr_two)

    fr_three = copy.deepcopy(base)
    fr_three["analysis_variant"] = (
        "FR three-component low-q candidate; norm-corrected nine-frame series"
    )
    window = lamellar_window(fr_three)
    original_broad = copy.deepcopy(window["peaks"][0])
    original_shoulder = copy.deepcopy(window["peaks"][1])
    window["peaks"] = [
        {
            "key": "q0295_broad",
            "label": "broad low-q component",
            "q_guess": 0.295,
            "q_min": 0.280,
            "q_max": 0.304,
            "fwhm_guess": 0.065,
            "fwhm_min": 0.025,
            "fwhm_max": 0.120,
            "anchor_temperature_C": 150,
            "provisional_assignment": True,
        },
        {
            **original_broad,
            "key": "q0307_narrow",
            "label": "0.307 narrow component",
            "q_guess": 0.307,
            "q_min": 0.304,
            "q_max": 0.316,
            "fwhm_guess": 0.015,
            "fwhm_min": 0.008,
            "fwhm_max": 0.040,
            "provisional_assignment": True,
        },
        original_shoulder,
    ]
    write_config("peakfit_config_FR_3component_candidate.json", fr_three)

    ip = copy.deepcopy(base)
    ip["analysis_variant"] = (
        "IP cut-specific simplified model; norm-corrected nine-frame series"
    )
    ip["windows"] = revised_ip_windows(ip)
    write_config("peakfit_config_IP_revised.json", ip)

    oop = copy.deepcopy(base)
    oop["analysis_variant"] = (
        "OOP one-broad-component low-q model; norm-corrected nine-frame series"
    )
    window = lamellar_window(oop)
    window["title"] = "OOP broad low-q lamellar candidate"
    window["peaks"] = [
        {
            "key": "q0306_broad_oop",
            "label": "OOP broad low-q component",
            "q_guess": 0.306,
            "q_min": 0.280,
            "q_max": 0.325,
            "fwhm_guess": 0.070,
            "fwhm_min": 0.025,
            "fwhm_max": 0.140,
            "anchor_temperature_C": 150,
            "provisional_assignment": True,
        }
    ]
    write_config("peakfit_config_OOP.json", oop)

    fr_revised = copy.deepcopy(fr_three)
    fr_revised["analysis_variant"] = (
        "FR revised low-degeneracy windows; norm-corrected nine-frame series"
    )
    fr_lamellar = lamellar_window(fr_revised)
    fr_lamellar["title"] = "FR descriptive three-component low-q model"
    fr_lamellar["peaks"][0]["anchor_temperature_C"] = 150
    fr_lamellar["peaks"][1]["anchor_temperature_C"] = 80
    fr_lamellar["peaks"][2]["anchor_temperature_C"] = 40
    fr_revised["windows"] = revised_windows(fr_revised)
    write_config("peakfit_config_FR_revised.json", fr_revised)

    oop_revised = copy.deepcopy(oop)
    oop_revised["analysis_variant"] = (
        "OOP revised low-degeneracy windows; norm-corrected nine-frame series"
    )
    oop_revised["windows"] = revised_windows(oop_revised)
    write_config("peakfit_config_OOP_revised.json", oop_revised)

    oop_quadratic = copy.deepcopy(oop_revised)
    oop_quadratic["analysis_variant"] = (
        "OOP low-q quadratic-background sensitivity check only"
    )
    sensitivity_window = copy.deepcopy(lamellar_window(oop_quadratic))
    sensitivity_window["background_order"] = 2
    sensitivity_window["title"] += " (quadratic-background sensitivity)"
    oop_quadratic["windows"] = [sensitivity_window]
    write_config(
        "peakfit_config_OOP_lowq_quadratic_sensitivity.json",
        oop_quadratic,
    )

    for cut, source in (
        ("FR", fr_revised),
        ("IP", ip),
        ("OOP", oop_revised),
    ):
        for two_components in (False, True):
            comparison = copy.deepcopy(source)
            comparison["analysis_variant"] = (
                f"{cut} high-temperature pi-pi component-count sensitivity"
            )
            comparison["windows"] = [
                high_temperature_pi_pi_window(two_components=two_components)
            ]
            suffix = "two" if two_components else "one"
            write_config(
                f"peakfit_config_{cut}_pipi_highT_{suffix}.json",
                comparison,
            )


if __name__ == "__main__":
    main()
