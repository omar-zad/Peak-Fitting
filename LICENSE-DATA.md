# Licence for the example data

The **software** in this repository (everything under `scripts/`, `tests/` and
the configuration files) is released under the MIT Licence in [LICENSE](LICENSE).

The **example data** is licensed separately.

## Scope

The following are data, not software:

- `example_data/drop40/` — processed one-dimensional FR, IP and OOP line cuts
  for the nine-frame Drop40 temperature series, plus the frame manifest.
- `example_data/secondary_scan/` — the frame manifest of scan 587178 (frame
  times, measured temperatures, phase labels and monitor counts; no
  intensities), used to exercise the configuration builders.
- `example_results/drop40_FR_checkpoint/` — the committed reference run derived
  from the Drop40 data.
- `example_results/drop40_IP_component_count/` — the committed component-count
  comparison derived from the Drop40 IP data.

## Terms

These files are released under the
[Creative Commons Attribution 4.0 International licence (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You may share and adapt them, including commercially, provided you give
appropriate credit, link to the licence, and indicate whether changes were
made. Cite this repository as described in [CITATION.cff](CITATION.cff).

## Provenance and limits

These are **processed, integrated line profiles**. The repository contains no
raw detector images, NeXus/HDF5 files, masks or `.poni` calibration files, and
no unrestricted source data.

The measurements were collected at a synchrotron beamline. Redistribution of
the underlying raw data remains subject to that facility's data policy and to
any applicable institutional agreement, independently of this licence. Anyone
wishing to obtain the raw data should contact the authors rather than assume
this licence covers it.

Sample labels in the secondary-sample configurations refer to collaborator
sample sets and are retained for provenance.
