# Changelog

## 1.1.0 (unreleased)

### Fixed
- **Eckart tunnelling.** The barrier parameters C and δ were wrong: C = (hν)² / [8(√V₁+√V₂)²], and δ used
  4V₁V₂ instead of B = (√V₁+√V₂)². The barrier behaved as almost transparent, and κ_Eckart(298 K) came out
  ≈ 5000 for the demo reaction. The code now uses C = (hν)²B/(16V₁V₂) (Johnston and Heicklen, 1962). P(E) is
  validated against a numerical solution of the 1D Schrödinger equation (`tests/test_eckart_reference.py`).
- **ORCA NEB support.** The CLI imported the ORCA parser but never called it, and the parser expected a
  format that ORCA does not write. The new parser reads ORCA 5/6 NEB, NEB-CI and NEB-TS outputs:
  - the `PATH SUMMARY` table (energies, distances, max|Fp| converted to eV/Å);
  - the climbing image and the NEB convergence banner;
  - the transition state from `PATH SUMMARY FOR NEB-TS`, with its convergence;
  - the final vibrational frequencies.

  `nebcert assess -i run.out` uses it automatically. When an optimised NEB-TS transition state is available,
  its barrier is used in preference to the spline maximum of the band.
- The climbing-image force warning now reports the force on the climbing image, not the band maximum.
- The demo output is labelled as synthetic data.

### Validation (`validation/orca_runs/`, ORCA 6.1.1, B3LYP-D3BJ/def2-SVP)
- **NH₃ inversion.** NEB-TS barrier 5.94 kcal/mol against 5.93 kcal/mol from the independently optimised
  planar TS; imaginary mode −830.15 cm⁻¹ against −830.27 cm⁻¹.
- **HCN → HNC.** Barrier 47.8 kcal/mol, HNC 13.6 kcal/mol above HCN, one imaginary mode (−1116 cm⁻¹).

### Not yet validated
- The Gaussian IRC parser has not been tested on real Gaussian output and is not used by the CLI.

## 1.0.0 (2026-08-31)

Initial release.
