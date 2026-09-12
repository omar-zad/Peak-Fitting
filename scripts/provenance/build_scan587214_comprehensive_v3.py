#!/usr/bin/env python3
"""Build comprehensive, regime-aware candidate configs for scan 587214 (Drop 50), version 3.

Version 3 revises the v2 candidate set after visual QC of the v2 selected-frame
fits (results/dimitar_drop50/scan_587214/comprehensive_frame_fits_v2):

* The low-q region is split into three regimes.  Frames 0-30 (early 46 C hold)
  show a very sharp line at 0.343 A^-1 with a shoulder near 0.357; by frame 50
  this has become a single peak near 0.350.  The fitter's new ``frame_min`` /
  ``frame_max`` window gate separates the early hold (frames 0-30) from the
  late hold (frames 31-52); the >= 65 C regime is unchanged.
* Lines that are only 1-2 q-steps wide (1.65 A^-1 at high temperature, the
  early 0.343 line, the IP pre-heating 0.35 line) are declared position_only
  so they are fitted freely instead of falling back to a borrowed shape.
* The high-temperature pi window separates the sharp conditional 1.85 line
  from the broad high-q shoulder near 1.87 that develops during the hold, and
  allows the main hump to be broader.
* FR near 0.49 gains the conditional 0.53 shoulder already used in OOP, so the
  asymmetric frame-99 profile can be fitted.
* IP near 0.58: the upper centre bound is lowered to 0.600 so the component
  cannot latch onto a two-point spike at 0.611-0.617 in frame 57.
* Windows are retitled where the v2 titles implied a temperature regime that
  the data do not support (the 0.60/0.67 features are already present at
  frame 50).
* OOP near 2.19: the second line at frame 50 sat exactly on the former lower
  centre bound (2.190) and was therefore refused by the gate; the bounds of the
  2.17/2.19 pair were separated at 2.181/2.182 (post-QC correction).
* FR: the 0.53 shoulder (frame 99) and the 1.96 line (frame 60) came to rest
  close to their lower centre bounds. The 1.96 range was widened to 1.935. For
  the shoulder a test with q_min 0.505 returned the identical frame-99 centre
  (0.5153, a free optimum) but let the shoulder steal the frame-85 peak flank,
  so the bound is set to 0.512 (post-QC correction).


Why this exists
---------------
The v1 selected-frame package (``expanded_locked_v1``) inherited six narrow
windows that had been tuned around frame 57.  Regions near 0.60, 0.67, 0.74,
1.50, 1.65, 1.75, 1.86, 2.03, 2.08, 2.17 and 2.27 A^-1 were never tested, and
the high-temperature broad pi-region hump near 1.74 A^-1 could not be fitted
because its centre was locked to 1.688-1.710 A^-1 inside a window that ended
at 1.735 A^-1.  These v2 configs restore scan-specific candidate coverage for
FR, IP and OOP.

Evidence base
-------------
Every window below was placed from the unsmoothed d5i-normalised profiles of
the eight selected frames (0, 25, 50, 57, 60, 71, 85, 99) using a
Savitzky-Golay local-maximum scan plus visual inspection.  ``evidence_frames``
records the selected frames in which a candidate was visible by eye; it is
informational only and does not influence the fit.

Modelling rules
---------------
* Pre-heating (<= 64.9 C, frames 0/25/50) and high-temperature (>= 65 C,
  frames 57-99) pi-region windows are separated because the sharp early lines
  are replaced by a broad hump during the 76 C hold.
* The high-temperature pi window carries both the sharp lines (still present
  at frames 57-71 and, in OOP, at frame 85) and a broad component, so the
  hump at frames 85/99 is fitted rather than pushed into the background.
* Main peaks are ordinary components so that their integrated areas remain
  reportable.  Lines that are only 1-3 q-steps wide are ``position_only``:
  their q0 may be reported but their width and area are suppressed.
* Nothing here is a crystallographic assignment, and no FWHM/CCL claim is
  intended.

Deliberately not fitted
-----------------------
* IP maximum near 1.046 A^-1: present with constant q in every frame including
  the amorphous frame 0 and absent in FR/OOP; treated as a fixed
  detector/substrate feature.
* A broad maximum drifting from ~1.17 to ~1.29 A^-1 over frames 0-32 in all
  cuts, then vanishing: a moving liquid/amorphous halo, not a lattice peak.
* OOP feature near 0.83 A^-1: alternates between a dip and a spike across
  frames, consistent with a detector-gap artefact.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "configs" / "in_situ" / "scan_587214" / "scan587214_FR.json"
OUTPUT = ROOT / "configs" / "in_situ" / "scan_587214" / "comprehensive_v3"
SELECTED_FRAMES = [0, 25, 50, 57, 60, 71, 85, 99]
FRAME_TEMPERATURES = {0: 45.9, 25: 46.5, 50: 46.6, 57: 77.0, 60: 76.1, 71: 76.1, 85: 76.3, 99: 76.1}
PROFILES = ["gaussian", "lorentzian", "pseudo_voigt"]
LOW_T_MAX = 64.9
HIGH_T_MIN = 65.0
EARLY_HOLD_LAST_FRAME = 30
LATE_HOLD_FIRST_FRAME = 31
LATE_HOLD_LAST_FRAME = 52
Q_STEP = 0.00622


def peak(
    key: str,
    label: str,
    guess: float,
    q_min: float,
    q_max: float,
    fwhm_guess: float,
    fwhm_min: float,
    fwhm_max: float,
    anchor_frame: int,
    *,
    evidence: list[int],
    position_only: bool = False,
    conditional: bool = False,
) -> dict[str, Any]:
    if not (q_min < guess < q_max):
        raise ValueError(f"{key}: guess outside bounds")
    if not (0 < fwhm_min <= fwhm_guess <= fwhm_max):
        raise ValueError(f"{key}: bad FWHM bounds")
    result: dict[str, Any] = {
        "key": key,
        "label": label,
        "q_guess": guess,
        "q_min": q_min,
        "q_max": q_max,
        "fwhm_guess": fwhm_guess,
        "fwhm_min": fwhm_min,
        "fwhm_max": fwhm_max,
        "anchor_frame": anchor_frame,
        "anchor_temperature_C": FRAME_TEMPERATURES[anchor_frame],
        "provisional_assignment": True,
        "evidence_frames": list(evidence),
        "conditional_candidate": bool(conditional),
    }
    if position_only:
        result["position_only"] = True
        result["suppress_fwhm_interpretation"] = True
    return result


def window(
    key: str,
    title: str,
    q_min: float,
    q_max: float,
    peaks: list[dict[str, Any]],
    *,
    temperature_min: float | None = None,
    temperature_max: float | None = None,
    frame_min: int | None = None,
    frame_max: int | None = None,
    note: str | None = None,
) -> dict[str, Any]:
    n_points = int(round((q_max - q_min) / Q_STEP)) + 1
    n_params = 2 + 3 * len(peaks)
    if n_points < n_params + 6:
        raise ValueError(f"{key}: only {n_points} points for {n_params} parameters")
    result: dict[str, Any] = {
        "key": key,
        "title": title,
        "q_min": q_min,
        "q_max": q_max,
        "background_order": 1,
        "peaks": peaks,
        "profiles_to_compare": list(PROFILES),
        "candidate_window_requires_visual_review": True,
    }
    if temperature_min is not None:
        result["temperature_min_C"] = temperature_min
    if temperature_max is not None:
        result["temperature_max_C"] = temperature_max
    if frame_min is not None:
        result["frame_min"] = frame_min
    if frame_max is not None:
        result["frame_max"] = frame_max
    if note:
        result["scientific_note"] = note
    return result


# ---------------------------------------------------------------------------
# FR
# ---------------------------------------------------------------------------
def fr_windows() -> list[dict[str, Any]]:
    return [
        window(
            "fr_low_q_early_hold", "FR low-q region, early 46 C hold (frames 0-30): broad component, sharp 0.343 line and 0.36 shoulder", 0.235, 0.420,
            [
                peak("fr_q0300_broad_early", "FR broad component near 0.30 (early hold)", 0.298, 0.272, 0.325, 0.055, 0.025, 0.130, 25, evidence=[25]),
                peak("fr_q0343_early", "FR very sharp line near 0.343 (early hold; position only)", 0.343, 0.334, 0.352, 0.012, 0.006, 0.030, 25, evidence=[25], position_only=True),
                peak("fr_q0358_early", "FR shoulder near 0.36 (early hold)", 0.358, 0.350, 0.375, 0.025, 0.014, 0.070, 25, evidence=[25]),
            ],
            temperature_max=LOW_T_MAX, frame_max=EARLY_HOLD_LAST_FRAME,
            note="Frames 22-30 show a one-to-two-pixel line at 0.343 with a shoulder near 0.357; both merge into one peak near 0.350 by frame 35. Frame 0 is amorphous.",
        ),
        window(
            "fr_low_q_late_hold", "FR low-q region, late 46 C hold (frames 31-52): broad component and 0.35 peak", 0.235, 0.420,
            [
                peak("fr_q0300_broad_late", "FR broad component near 0.30 (late hold)", 0.300, 0.272, 0.325, 0.055, 0.025, 0.130, 50, evidence=[50]),
                peak("fr_q0350_late", "FR sharper component near 0.35 (late hold)", 0.350, 0.332, 0.370, 0.022, 0.009, 0.060, 50, evidence=[50]),
            ],
            temperature_max=LOW_T_MAX, frame_min=LATE_HOLD_FIRST_FRAME, frame_max=LATE_HOLD_LAST_FRAME,
        ),
        window(
            "fr_low_q_highT", "FR low-q region, 76 C hold: broad component and 0.35 peak", 0.235, 0.420,
            [
                peak("fr_q0300_broad", "FR broad component near 0.30", 0.300, 0.272, 0.325, 0.055, 0.025, 0.130, 71, evidence=[57, 60, 71, 85, 99]),
                peak("fr_q0350", "FR sharper component near 0.35", 0.350, 0.332, 0.370, 0.022, 0.009, 0.060, 60, evidence=[57, 60, 71, 85, 99]),
            ],
            temperature_min=HIGH_T_MIN,
        ),
        window(
            "fr_q0495_q0530", "FR component near 0.49 and conditional shoulder near 0.53", 0.435, 0.565,
            [
                peak("fr_q0495", "FR component near 0.49", 0.494, 0.462, 0.522, 0.030, 0.010, 0.100, 60, evidence=[25, 50, 57, 60, 71, 85, 99]),
                peak("fr_q0530", "FR conditional shoulder near 0.52-0.53", 0.525, 0.512, 0.550, 0.030, 0.009, 0.080, 99, evidence=[85, 99], conditional=True),
            ],
            note="At frames 85-99 the 0.49 maximum moves to ~0.486 and a shoulder near 0.52-0.53 appears; a single symmetric profile no longer fits.",
        ),
        window(
            "fr_q0600_q0665", "FR features near 0.60 and 0.67 (present from the late 46 C hold onward)", 0.560, 0.700,
            [
                peak("fr_q0600", "FR component near 0.60", 0.600, 0.575, 0.628, 0.030, 0.009, 0.120, 99, evidence=[50, 57, 60, 71, 85, 99]),
                peak("fr_q0665", "FR component near 0.67", 0.665, 0.640, 0.688, 0.025, 0.009, 0.080, 85, evidence=[50, 71, 85, 99]),
            ],
            note="Both features are already weakly visible at frame 50 (end of the 46 C hold); their onset lies between frames 25 and 50.",
        ),
        window(
            "fr_q0735", "FR sharp line 0.725-0.738 (0.738 at 46 C, 0.725-0.729 at 76 C)", 0.695, 0.790,
            [peak("fr_q0735", "FR sharp component near 0.74 (position only)", 0.735, 0.715, 0.760, 0.018, 0.008, 0.060, 50, evidence=[50, 57, 60], position_only=True)],
        ),
        window(
            "fr_q0937", "FR sharp feature near 0.94", 0.880, 1.000,
            [peak("fr_q0937", "FR sharp component near 0.94 (position only)", 0.937, 0.918, 0.958, 0.016, 0.008, 0.060, 50, evidence=[25, 50, 57, 60, 71], position_only=True)],
        ),
        window(
            "fr_q1460_q1500_q1536", "FR lines near 1.46, 1.50 and 1.54", 1.405, 1.600,
            [
                peak("fr_q1460", "FR conditional line near 1.46 (position only)", 1.460, 1.440, 1.480, 0.018, 0.008, 0.060, 50, evidence=[50], position_only=True, conditional=True),
                peak("fr_q1500", "FR line near 1.50 (position only)", 1.500, 1.482, 1.516, 0.020, 0.008, 0.070, 50, evidence=[50, 57, 60], position_only=True),
                peak("fr_q1536", "FR sharp line near 1.54 (position only)", 1.536, 1.522, 1.552, 0.014, 0.008, 0.050, 50, evidence=[25, 50, 57, 60, 71], position_only=True),
            ],
        ),
        window(
            "fr_pi_lowT", "FR pre-heating pi-region lines", 1.610, 1.920,
            [
                peak("fr_q1660_lowT", "FR pre-heating line near 1.66 (position only)", 1.660, 1.640, 1.680, 0.016, 0.008, 0.060, 50, evidence=[50], position_only=True),
                peak("fr_q1708_lowT", "FR pre-heating line near 1.71 (position only)", 1.708, 1.690, 1.725, 0.016, 0.008, 0.060, 50, evidence=[25, 50], position_only=True),
                peak("fr_q1752_lowT", "FR pre-heating line near 1.75 (position only)", 1.752, 1.735, 1.772, 0.022, 0.008, 0.080, 50, evidence=[25, 50], position_only=True),
                peak("fr_q1858_lowT", "FR conditional pre-heating line near 1.86 (position only)", 1.858, 1.835, 1.882, 0.025, 0.008, 0.090, 50, evidence=[25, 50], position_only=True, conditional=True),
            ],
            temperature_max=LOW_T_MAX,
            note="Frame 0 is amorphous in this region; components are expected to be not_detected there.",
        ),
        window(
            "fr_pi_highT", "FR high-temperature pi-region: sharp lines plus emerging broad hump", 1.600, 1.940,
            [
                peak("fr_q1650_highT", "FR high-temperature sharp line near 1.65 (position only)", 1.650, 1.630, 1.672, 0.014, 0.008, 0.060, 60, evidence=[57, 60, 71], position_only=True),
                peak("fr_q1700_highT", "FR high-temperature sharp line near 1.70 (position only)", 1.699, 1.684, 1.716, 0.014, 0.008, 0.050, 60, evidence=[57, 60, 71, 85], position_only=True),
                peak("fr_q1738_highT", "FR high-temperature sharp shoulder near 1.74 (position only)", 1.738, 1.724, 1.756, 0.022, 0.008, 0.050, 60, evidence=[57, 60, 71], position_only=True),
                peak("fr_q1750_broad_highT", "FR broad high-temperature hump near 1.73-1.75", 1.745, 1.700, 1.800, 0.130, 0.060, 0.350, 99, evidence=[57, 60, 71, 85, 99]),
                peak("fr_q1856_highT", "FR conditional sharp line near 1.86 (position only)", 1.856, 1.840, 1.872, 0.012, 0.008, 0.035, 60, evidence=[57, 60], position_only=True, conditional=True),
                peak("fr_q1870_shoulder_highT", "FR broad high-q shoulder near 1.85-1.88", 1.865, 1.825, 1.905, 0.080, 0.040, 0.160, 71, evidence=[57, 60, 71, 85, 99]),
            ],
            temperature_min=HIGH_T_MIN,
            note="The broad hump and the broad high-q shoulder are background-sensitive by construction; their centres are descriptive positions, not lattice spacings.",
        ),
        window(
            "fr_q1965_q2008_q2030", "FR lines near 1.94-1.96, 2.01 and 2.03", 1.915, 2.060,
            [
                peak("fr_q1965", "FR conditional line near 1.94-1.96 (position only)", 1.958, 1.935, 1.980, 0.018, 0.008, 0.060, 50, evidence=[50, 57, 60], position_only=True, conditional=True),
                peak("fr_q2008", "FR line near 2.01 (position only)", 2.008, 1.994, 2.020, 0.018, 0.008, 0.060, 50, evidence=[25, 50, 57, 60, 71], position_only=True),
                peak("fr_q2030", "FR conditional line near 2.03 (position only)", 2.030, 2.021, 2.042, 0.012, 0.008, 0.040, 50, evidence=[50], position_only=True, conditional=True),
            ],
        ),
        window(
            "fr_q2085_q2172", "FR lines near 2.08 and 2.17", 2.050, 2.215,
            [
                peak("fr_q2085", "FR conditional line near 2.08 (position only)", 2.085, 2.062, 2.105, 0.020, 0.008, 0.070, 50, evidence=[50], position_only=True, conditional=True),
                peak("fr_q2172", "FR line near 2.17 (position only)", 2.172, 2.152, 2.190, 0.020, 0.008, 0.070, 50, evidence=[50, 57, 60], position_only=True),
            ],
        ),
        window(
            "fr_q2270", "FR conditional line near 2.25-2.28", 2.220, 2.330,
            [peak("fr_q2270", "FR conditional line near 2.27 (position only)", 2.268, 2.240, 2.295, 0.020, 0.008, 0.080, 50, evidence=[25, 50], position_only=True, conditional=True)],
        ),
    ]


# ---------------------------------------------------------------------------
# IP
# ---------------------------------------------------------------------------
def ip_windows() -> list[dict[str, Any]]:
    return [
        window(
            "ip_low_q_lowT", "IP pre-heating low-q components", 0.225, 0.415,
            [
                peak("ip_q0295_broad_lowT", "IP pre-heating broad component near 0.30-0.33", 0.300, 0.258, 0.340, 0.070, 0.030, 0.140, 50, evidence=[50]),
                peak("ip_q0348_lowT", "IP pre-heating sharp line near 0.35 (position only)", 0.348, 0.333, 0.364, 0.014, 0.006, 0.050, 50, evidence=[25, 50], position_only=True),
            ],
            temperature_max=LOW_T_MAX,
        ),
        window(
            "ip_low_q_highT", "IP high-temperature low-q components", 0.225, 0.415,
            [
                peak("ip_q0287_broad_highT", "IP high-temperature broad component near 0.29", 0.287, 0.258, 0.315, 0.070, 0.030, 0.140, 99, evidence=[71, 85, 99]),
                peak("ip_q0347_highT", "IP high-temperature component near 0.35", 0.347, 0.333, 0.364, 0.025, 0.009, 0.070, 85, evidence=[57, 60, 71, 85, 99]),
                peak("ip_q0387_highT", "IP late sharp line near 0.39 (position only)", 0.387, 0.376, 0.398, 0.012, 0.008, 0.040, 99, evidence=[71, 85, 99], position_only=True),
            ],
            temperature_min=HIGH_T_MIN,
        ),
        window(
            "ip_q0474_q0505", "IP components near 0.47 and 0.50", 0.430, 0.545,
            [
                peak("ip_q0474", "IP line near 0.47 (position only)", 0.474, 0.458, 0.487, 0.016, 0.008, 0.060, 57, evidence=[57], position_only=True),
                peak("ip_q0505", "IP component near 0.50", 0.505, 0.488, 0.522, 0.025, 0.009, 0.080, 99, evidence=[50, 57, 60, 71, 85, 99]),
            ],
        ),
        window(
            "ip_q0585", "IP late component near 0.58", 0.545, 0.645,
            [peak("ip_q0585", "IP component near 0.58", 0.582, 0.555, 0.600, 0.030, 0.009, 0.090, 99, evidence=[71, 85, 99])],
        ),
        window(
            "ip_q0750", "IP conditional line near 0.75", 0.700, 0.805,
            [peak("ip_q0750", "IP conditional line near 0.75 (position only)", 0.750, 0.732, 0.768, 0.016, 0.008, 0.060, 50, evidence=[50, 57], position_only=True, conditional=True)],
        ),
        window(
            "ip_q1502_q1543", "IP lines near 1.50 and 1.54", 1.455, 1.585,
            [
                peak("ip_q1502", "IP component near 1.50", 1.502, 1.485, 1.520, 0.025, 0.009, 0.080, 50, evidence=[50, 57, 60]),
                peak("ip_q1543", "IP sharp line near 1.54 (position only)", 1.543, 1.528, 1.558, 0.014, 0.008, 0.050, 25, evidence=[25, 50, 57, 60, 71], position_only=True),
            ],
        ),
        window(
            "ip_pi_lowT", "IP pre-heating pi-region components", 1.620, 1.930,
            [
                peak("ip_q1668_lowT", "IP conditional pre-heating line near 1.67 (position only)", 1.668, 1.648, 1.688, 0.018, 0.008, 0.070, 50, evidence=[50], position_only=True, conditional=True),
                peak("ip_q1755_lowT", "IP pre-heating component near 1.75", 1.755, 1.730, 1.780, 0.030, 0.009, 0.110, 50, evidence=[25, 50]),
                peak("ip_q1875_lowT", "IP conditional pre-heating line near 1.88 (position only)", 1.875, 1.855, 1.895, 0.020, 0.008, 0.080, 50, evidence=[50], position_only=True, conditional=True),
            ],
            temperature_max=LOW_T_MAX,
        ),
        window(
            "ip_pi_highT", "IP high-temperature pi-region: broad hump, broad high-q shoulder and conditional 1.66 line", 1.580, 1.940,
            [
                peak("ip_q1660_highT", "IP conditional high-temperature line near 1.66 (position only)", 1.660, 1.645, 1.678, 0.014, 0.008, 0.050, 57, evidence=[57, 60], position_only=True, conditional=True),
                peak("ip_q1750_broad_highT", "IP broad high-temperature hump near 1.75", 1.750, 1.705, 1.800, 0.100, 0.030, 0.300, 99, evidence=[57, 60, 71, 85, 99]),
                peak("ip_q1875_shoulder_highT", "IP broad high-q shoulder near 1.87", 1.875, 1.830, 1.910, 0.080, 0.040, 0.180, 99, evidence=[57, 60, 71, 85, 99]),
            ],
            temperature_min=HIGH_T_MIN,
        ),
        window(
            "ip_q2016", "IP sharp line near 2.02", 1.955, 2.055,
            [peak("ip_q2016", "IP sharp line near 2.02 (position only)", 2.016, 1.998, 2.030, 0.014, 0.008, 0.050, 50, evidence=[25, 50, 57, 60, 71], position_only=True)],
        ),
        window(
            "ip_q2080_q2115", "IP lines near 2.08 and 2.12", 2.045, 2.160,
            [
                peak("ip_q2080", "IP line near 2.08 (position only)", 2.080, 2.062, 2.096, 0.016, 0.008, 0.060, 50, evidence=[50, 57, 60], position_only=True),
                peak("ip_q2115", "IP line near 2.12 (position only)", 2.115, 2.100, 2.132, 0.016, 0.008, 0.060, 50, evidence=[50, 57, 60], position_only=True),
            ],
        ),
        window(
            "ip_q2254", "IP conditional line near 2.25", 2.205, 2.310,
            [peak("ip_q2254", "IP conditional line near 2.25 (position only)", 2.254, 2.238, 2.272, 0.014, 0.008, 0.060, 25, evidence=[25, 50], position_only=True, conditional=True)],
        ),
    ]


# ---------------------------------------------------------------------------
# OOP
# ---------------------------------------------------------------------------
def oop_windows() -> list[dict[str, Any]]:
    return [
        window(
            "oop_low_q_early_hold", "OOP low-q region, early 46 C hold (frames 0-30): broad component, sharp 0.343 line and 0.36 shoulder", 0.240, 0.420,
            [
                peak("oop_q0300_broad_early", "OOP broad component near 0.30 (early hold)", 0.298, 0.272, 0.328, 0.060, 0.025, 0.140, 25, evidence=[25]),
                peak("oop_q0343_early", "OOP very sharp line near 0.343 (early hold; position only)", 0.343, 0.334, 0.352, 0.010, 0.006, 0.030, 25, evidence=[25], position_only=True),
                peak("oop_q0360_early", "OOP shoulder near 0.36 (early hold)", 0.360, 0.350, 0.378, 0.025, 0.014, 0.070, 25, evidence=[25]),
            ],
            temperature_max=LOW_T_MAX, frame_max=EARLY_HOLD_LAST_FRAME,
            note="Frames 22-30 show a one-pixel line at 0.343 with a shoulder near 0.36; both merge into one peak near 0.353 by frame 35. Frame 0 is amorphous.",
        ),
        window(
            "oop_low_q_late_hold", "OOP low-q region, late 46 C hold (frames 31-52): broad component and 0.35 peak", 0.240, 0.420,
            [
                peak("oop_q0305_broad_late", "OOP broad component near 0.31 (late hold)", 0.308, 0.278, 0.328, 0.060, 0.025, 0.140, 50, evidence=[50]),
                peak("oop_q0355_late", "OOP sharper component near 0.35 (late hold)", 0.354, 0.336, 0.372, 0.022, 0.009, 0.065, 50, evidence=[50]),
            ],
            temperature_max=LOW_T_MAX, frame_min=LATE_HOLD_FIRST_FRAME, frame_max=LATE_HOLD_LAST_FRAME,
        ),
        window(
            "oop_low_q_highT", "OOP low-q region, 76 C hold: broad component and 0.35-0.36 peak", 0.240, 0.420,
            [
                peak("oop_q0305_broad", "OOP broad component near 0.31", 0.305, 0.278, 0.328, 0.060, 0.025, 0.140, 71, evidence=[57, 60, 71, 85, 99]),
                peak("oop_q0355", "OOP sharper component near 0.35-0.36", 0.354, 0.336, 0.372, 0.022, 0.009, 0.065, 60, evidence=[57, 60, 71, 85, 99]),
            ],
            temperature_min=HIGH_T_MIN,
        ),
        window(
            "oop_q0492_q0535", "OOP component near 0.49 and conditional shoulder near 0.54", 0.430, 0.565,
            [
                peak("oop_q0492", "OOP component near 0.49", 0.493, 0.470, 0.512, 0.028, 0.009, 0.100, 60, evidence=[25, 50, 57, 60, 71, 85, 99]),
                peak("oop_q0535", "OOP conditional shoulder near 0.54", 0.535, 0.518, 0.552, 0.030, 0.009, 0.080, 99, evidence=[71, 85, 99], conditional=True),
            ],
            note="At frames 71-99 the main peak becomes asymmetric with a sub-maximum near 0.480 two pixels below the 0.493 maximum; this is not resolvable on the present q grid and the fitted centre is the intensity-weighted maximum.",
        ),
        window(
            "oop_q0604_q0672", "OOP components near 0.60 and 0.67", 0.560, 0.705,
            [
                peak("oop_q0604", "OOP component near 0.60", 0.603, 0.578, 0.622, 0.030, 0.009, 0.090, 99, evidence=[25, 50, 57, 60, 71, 85, 99]),
                peak("oop_q0672", "OOP component near 0.67", 0.672, 0.650, 0.692, 0.028, 0.009, 0.085, 99, evidence=[60, 71, 85, 99]),
            ],
        ),
        window(
            "oop_q0736", "OOP sharp feature near 0.74", 0.695, 0.790,
            [peak("oop_q0736", "OOP sharp component near 0.74 (position only)", 0.736, 0.716, 0.756, 0.018, 0.008, 0.060, 50, evidence=[50, 57, 60, 71], position_only=True)],
        ),
        window(
            "oop_q0938", "OOP sharp feature near 0.94", 0.880, 1.000,
            [peak("oop_q0938", "OOP sharp component near 0.94 (position only)", 0.937, 0.918, 0.958, 0.016, 0.008, 0.060, 50, evidence=[25, 50, 57, 60, 71, 85, 99], position_only=True)],
        ),
        window(
            "oop_q1464_q1495_q1538", "OOP lines near 1.46, 1.50 and 1.54", 1.405, 1.595,
            [
                peak("oop_q1464", "OOP conditional line near 1.46 (position only)", 1.464, 1.445, 1.482, 0.020, 0.008, 0.070, 50, evidence=[50], position_only=True, conditional=True),
                peak("oop_q1495", "OOP conditional line near 1.50 (position only)", 1.495, 1.483, 1.510, 0.016, 0.008, 0.060, 25, evidence=[25], position_only=True, conditional=True),
                peak("oop_q1538", "OOP sharp line near 1.54 (position only)", 1.537, 1.522, 1.553, 0.014, 0.008, 0.050, 50, evidence=[50], position_only=True),
            ],
        ),
        window(
            "oop_pi_lowT", "OOP pre-heating pi-region lines", 1.610, 1.920,
            [
                peak("oop_q1657_lowT", "OOP pre-heating line near 1.66 (position only)", 1.657, 1.638, 1.676, 0.016, 0.008, 0.060, 50, evidence=[50], position_only=True),
                peak("oop_q1706_lowT", "OOP pre-heating sharp line near 1.71 (position only)", 1.706, 1.690, 1.722, 0.014, 0.008, 0.050, 50, evidence=[25, 50], position_only=True),
                peak("oop_q1745_lowT", "OOP pre-heating line near 1.74 (position only)", 1.745, 1.730, 1.768, 0.020, 0.008, 0.080, 50, evidence=[25, 50], position_only=True),
                peak("oop_q1856_lowT", "OOP pre-heating line near 1.86 (position only)", 1.856, 1.835, 1.878, 0.022, 0.008, 0.085, 50, evidence=[25, 50], position_only=True),
            ],
            temperature_max=LOW_T_MAX,
            note="Frame 0 is amorphous in this region; components are expected to be not_detected there.",
        ),
        window(
            "oop_pi_highT", "OOP high-temperature pi-region: sharp lines plus emerging broad hump", 1.600, 1.940,
            [
                peak("oop_q1650_highT", "OOP high-temperature sharp line near 1.65 (position only)", 1.650, 1.632, 1.670, 0.014, 0.008, 0.060, 60, evidence=[57, 60, 71], position_only=True),
                peak("oop_q1700_highT", "OOP high-temperature sharp line near 1.70 (position only)", 1.700, 1.685, 1.715, 0.014, 0.008, 0.050, 60, evidence=[57, 60, 71, 85], position_only=True),
                peak("oop_q1737_highT", "OOP high-temperature sharp shoulder near 1.74 (position only)", 1.737, 1.724, 1.755, 0.022, 0.008, 0.050, 60, evidence=[57, 60, 71, 85], position_only=True),
                peak("oop_q1750_broad_highT", "OOP broad high-temperature hump near 1.73-1.75", 1.745, 1.700, 1.800, 0.130, 0.060, 0.350, 99, evidence=[57, 60, 71, 85, 99]),
                peak("oop_q1856_highT", "OOP conditional sharp line near 1.86 (position only)", 1.856, 1.840, 1.872, 0.012, 0.008, 0.035, 60, evidence=[25, 50, 57, 60, 71], position_only=True, conditional=True),
                peak("oop_q1870_shoulder_highT", "OOP broad high-q shoulder near 1.85-1.88", 1.865, 1.825, 1.905, 0.080, 0.040, 0.160, 85, evidence=[71, 85, 99]),
            ],
            temperature_min=HIGH_T_MIN,
            note="The broad hump and the broad high-q shoulder are background-sensitive by construction; their centres are descriptive positions, not lattice spacings.",
        ),
        window(
            "oop_q2012_q2032", "OOP lines near 2.01 and 2.03", 1.945, 2.055,
            [
                peak("oop_q2012", "OOP line near 2.01 (position only)", 2.012, 1.996, 2.024, 0.016, 0.008, 0.060, 50, evidence=[25, 50], position_only=True),
                peak("oop_q2032", "OOP conditional line near 2.03 (position only)", 2.032, 2.024, 2.045, 0.012, 0.008, 0.040, 50, evidence=[50], position_only=True, conditional=True),
            ],
        ),
        window(
            "oop_q2168_q2195", "OOP lines near 2.17 and 2.19", 2.090, 2.235,
            [
                peak("oop_q2168", "OOP line near 2.17 (position only)", 2.166, 2.150, 2.181, 0.022, 0.008, 0.080, 50, evidence=[50, 57, 60], position_only=True),
                peak("oop_q2195", "OOP second line near 2.19 (position only)", 2.193, 2.182, 2.216, 0.018, 0.008, 0.060, 50, evidence=[50, 57], position_only=True),
            ],
        ),
    ]


def build(cut: str, windows: list[dict[str, Any]]) -> dict[str, Any]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    config: dict[str, Any] = {
        "frame_temperatures_C": source["frame_temperatures_C"],
        "profiles_to_compare": list(PROFILES),
        "fit_settings": dict(source["fit_settings"]),
        "windows": windows,
        "analysis_scope": (
            "Scan-specific comprehensive selected-frame screening for peak presence "
            "and q0/d across FR, IP and OOP. Position-only lines cannot supply FWHM, "
            "CCL or integrated-area claims."
        ),
        "representative_frames": list(SELECTED_FRAMES),
        "series_metadata": {
            "scan": 587214,
            "sample": "Dimitar_drop50",
            "cut": cut,
            "candidate_version": 3,
            "selected_frames": list(SELECTED_FRAMES),
            "source_of_windows": (
                "Unsmoothed scan-587214 selected-frame profiles (Savitzky-Golay "
                "maximum scan plus visual inspection); not transferred from another sample."
            ),
            "superseded_configs": [
                str(ROOT / "configs" / "in_situ" / "scan_587214" / "expanded_locked_v1" / f"scan587214_{cut}_locked_expanded.json"),
                str(ROOT / "configs" / "in_situ" / "scan_587214" / "comprehensive_v2" / f"scan587214_{cut}_comprehensive_v3.json"),
            ],
            "v3_revision_basis": "Visual QC of the v2 selected-frame fits (per-window overlays and residuals) plus a numerical audit of gate failures on visible features.",
            "regime_rule": (
                "Pre-heating windows use temperature <= 64.9 C (frames 0, 25, 50); "
                "high-temperature windows use >= 65 C (frames 57-99). The low-q "
                "pre-heating region is further split by acquisition frame (early hold "
                "frames 0-30 vs late hold frames 31-52) because the pattern changes "
                "with time at constant temperature. Components are still accepted "
                "separately in each frame by leave-one-component-out delta BIC and "
                "area SNR."
            ),
            "width_policy": (
                "Ordinary components keep width and area; position_only lines may "
                "support q0/presence while their width and area are suppressed."
            ),
            "conditional_candidates": (
                "conditional_candidate=true marks weak or single-frame candidates that "
                "are tested but expected to fail the gate in most frames."
            ),
            "deliberately_not_fitted": [
                "IP fixed maximum near 1.046 A^-1 (present in amorphous frame 0; detector/substrate feature).",
                "Drifting halo 1.17-1.29 A^-1 over frames 0-32 in all cuts (liquid/amorphous, not a lattice peak).",
                "OOP dip/spike near 0.83 A^-1 (detector-gap artefact).",
            ],
        },
    }
    config["fit_settings"]["multistarts_model_selection"] = 5
    config["fit_settings"]["multistarts_final"] = 8
    return config


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for cut, windows in {"FR": fr_windows(), "IP": ip_windows(), "OOP": oop_windows()}.items():
        output = OUTPUT / f"scan587214_{cut}_comprehensive_v3.json"
        output.write_text(json.dumps(build(cut, windows), indent=2) + "\n", encoding="utf-8")
        n_peaks = sum(len(w["peaks"]) for w in windows)
        print(f"{output}  windows={len(windows)} components={n_peaks}")


if __name__ == "__main__":
    main()
