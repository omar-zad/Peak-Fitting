# Pre-publication checklist

The repository is private. The steps below still require human decisions and,
in two cases, other people's approval, before it is made public or cited from
the dissertation. Nothing in the code makes the repository public.

## Blocking, before the repository goes public

- [ ] **Fill in the copyright holder.** `LICENSE` and `CITATION.cff` contain
      placeholders. Use your full legal name. Check whether your institution
      requires itself to be named as a copyright holder for work produced
      during a funded studentship; many do.
- [ ] **Confirm the data policy permits publishing the example data.**
      `example_data/` and `example_results/` contain processed data derived from
      beamtime. Confirm with your supervisor and the facility's data manager
      that publishing processed 1D profiles and a frame manifest is permitted,
      and whether an embargo period applies. This is the one item that cannot
      be resolved from inside the repository.
- [ ] **Confirm collaborator consent for sample labelling.** Sample labels such
      as `Dimitar_bc` appear in the secondary-sample configurations and
      documentation. Confirm this is acceptable to the person named, or rename
      the sample slugs.
- [ ] **Supervisor sign-off** on making the repository public.

## Recommended before citing it in the dissertation

- [ ] **Confirm CI is green** on GitHub. `.github/workflows/ci.yml` runs the
      full test suite, including the slow Drop40 chain, on Python 3.10, 3.12
      and 3.13.
- [ ] **Tag the release** used by the dissertation, for example `v1.0.0`, so
      the cited state is unambiguous even after later changes.
- [ ] **Mint a DOI.** Enable the GitHub-Zenodo integration, then publish a
      GitHub release from the tag. Zenodo mints a DOI that survives renaming,
      transfer or deletion of the GitHub repository. Cite the DOI, not the bare
      GitHub URL; a dissertation is permanent and a URL is not.
- [ ] **Add the DOI** to `CITATION.cff`.

## What the tests verify

`tests/test_reference_checkpoint.py` re-runs the committed Drop40 FR checkpoint
and asserts unchanged profile-family selection, fitted centres, widths, areas,
BIC, detection status and reporting-gate membership. It was validated on
Linux under two independent environments and on macOS:

| Environment | numpy | pandas | scipy | Result vs committed reference |
|---|---|---|---|---|
| Python 3.13, Linux | 2.4.6 | 3.0.5 | 1.17.1 | identical selection; max abs q difference 1.1e-9 |
| Python 3.10, Linux | 2.2.6 | 2.3.3 | 1.15.3 | identical selection; max abs q difference 1.1e-9 |
| Python 3.13, macOS arm64 | 2.5.3 | 3.0.5 | 1.18.1 | identical selection and gates |

The committed reference itself was produced on macOS (arm64, Python 3.13,
numpy 2.5.1).

Three windows in the checkpoint are decided by a BIC margin below one unit
(`pi_pi_sharp_pair`, `q0754`, `q1999_single`). Observed cross-platform BIC
variation is about 7e-11, so exact assertions on family selection are safe by
roughly nine orders of magnitude. If the test ever fails on one of those three
windows, treat it as a genuine near-tie flip to investigate, not as noise.

`tests/test_config_builders.py` asserts that `build_peakfit_cut_configs.py`
regenerates every file under `configs/drop40/` exactly, and that the two other
builders run on the committed example manifest. `tests/test_drop40_pipeline.py`
(slow) runs the whole documented Drop40 chain and asserts the three BIC values
quoted in the README, which it also compares with the committed comparison
table.
