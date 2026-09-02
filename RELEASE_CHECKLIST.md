# Pre-publication checklist

This branch prepares the repository for public release and citation from a
dissertation. **Nothing here makes the repository public.** The steps below
still require human decisions and, in two cases, other people's approval.

## Blocking — must be done before the repository goes public

- [ ] **Fill in the copyright holder.** `LICENSE` and `CITATION.cff` contain
      placeholders. Use your full legal name. Check whether your institution
      requires itself to be named as a copyright holder for work produced
      during a funded studentship — many do.
- [ ] **Confirm the data policy permits publishing the example line cuts.**
      `example_data/` and `example_results/` contain processed data derived from
      beamtime. Confirm with your supervisor and the facility's data manager
      that publishing processed 1D profiles is permitted, and whether an embargo
      period applies. This is the one item that cannot be resolved from inside
      the repository.
- [ ] **Confirm collaborator consent for sample labelling.** Sample labels such
      as `Dimitar_bc` appear in 25+ configuration files. Confirm this is
      acceptable to the person named, or rename the sample slugs.
- [ ] **Supervisor sign-off** on making the repository public.

## Recommended before citing it in the dissertation

- [ ] **Mint a DOI.** Enable the GitHub–Zenodo integration, then publish a
      release. Zenodo mints a DOI that survives renaming, transfer or deletion
      of the GitHub repository. Cite the DOI, not the bare GitHub URL — a
      dissertation is permanent and a URL is not.
- [ ] **Add the DOI** to `CITATION.cff` and to the README badge line.
- [ ] **Confirm CI passes** on the first push (`.github/workflows/ci.yml` runs
      the reference-run regression test on Python 3.10, 3.12 and 3.13).
- [ ] **Tag the release** used by the dissertation, e.g. `v1.0.0`, so the cited
      state is unambiguous even after later changes.

## Notes on what was verified

The reference-run regression test in `tests/` was validated by re-running the
committed Drop40 FR checkpoint on Linux under two independent environments:

| Environment | numpy | pandas | scipy | Result vs committed reference |
|---|---|---|---|---|
| Python 3.13 | 2.4.6 | 3.0.5 | 1.17.1 | identical selection; max abs q difference 1.1e-9 |
| Python 3.10 | 2.2.6 | 2.3.3 | 1.15.3 | identical selection; max abs q difference 1.1e-9 |

The committed reference itself was produced on macOS (arm64, Python 3.13,
numpy 2.5.1). Profile-family selection, detection status and the set of
gate-passing reported positions were identical in every case.

Three windows in the checkpoint are decided by a BIC margin below one unit
(`pi_pi_sharp_pair`, `q0754`, `q1999_single`). Observed cross-platform BIC
variation is about 7e-11, so exact assertions on family selection are safe by
roughly nine orders of magnitude. If the test ever fails on one of those three
windows, treat it as a genuine near-tie flip to investigate, not as noise.
