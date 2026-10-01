# Changelog

## 1.2.0 (unreleased)

Version 1.1.0 was committed to the repository but never released; 1.2.0 contains its changes (listed under
1.1.0 below) and the following. Version 1.0.0 is the only earlier release.

### Fixed
- **`--irc` reported a failed IRC whenever it was omitted.** The flag was `store_true`, so leaving it out passed
  `irc_confirmed=False` and every transition state got a WARNING saying that the IRC did not connect the intended
  minima. Omitting the flag now means that no IRC was reported (None).
- **Convergence was read but never checked.** A band that ORCA did not report as converged, or a NEB-TS
  transition-state optimisation that did not converge, now gives FAIL. 1.0.0 and 1.1.0 passed them.
- **Failed ORCA jobs.** ORCA prints the `PATH SUMMARY` table once before the band is optimised, with E = 0 for
  the images; a job whose image calculations fail leaves only this table. 1.0.0 and 1.1.0 read it as a
  71,706 kcal/mol barrier and passed it; the table is now rejected with an explanation.
- **Identical end points.** When pre-optimisation moves one end point into the other, the band has zero length.
  It is now a FAIL; distances printed as 0.000 for every image no longer crash the spline (1.0.0 and 1.1.0 raised
  `ValueError`).
- **Barrier source.** The optimised NEB-TS energy is the barrier of the band check as well as of the kinetics;
  without one, the energy of a converged climbing image is used instead of the spline maximum, which overshot
  the saddle by 1.3 kcal/mol on a 3-image band. The spline barrier is kept for comparison
  (`e_forward_barrier_spline_ev`, `barrier_source`).
- **Band force check after a converged TS optimisation.** NEB-TS converges the band only loosely by design;
  the 0.05 eV/Å check on the image nearest the barrier is now skipped when the TS optimisation converged (it
  warned on four of eight converged NEB-TS runs).
- **Intermediate minima** are reported only when deeper than 0.5 kcal/mol below both bounding maxima
  (`min_well_depth_kcal_mol`); any dip between neighbouring images used to give a warning.
- **Imaginary frequencies.** Modes between −50 and −1 cm⁻¹ are counted separately and give a WARNING; 1.0.0 and
  1.1.0 counted them as a second imaginary mode and reported a higher-order saddle point (FAIL).
- **Eyring rate from an electronic barrier.** The expression was labelled a TST rate constant whatever the barrier.
  `calculate_eyring_tst_rates(..., barrier_type="electronic")` (the default) now reports the number with the
  status NOT_APPLICABLE and says that zero-point, thermal and entropic terms are missing; it does not set the
  overall status. `barrier_type="gibbs"` and the CLI option `--gibbs-barrier-kcal` give the Eyring rate
  constant.
- **Tunnelling regime.** A WARNING is given when 298.15 K is below the crossover temperature
  T_c = hc|ν|/(2πk_B), where the Wigner factor is invalid and the one-dimensional Eckart factor is only an
  estimate; a non-positive barrier gives FAIL.
- **Wording.** "FULLY VALIDATED (PUBLICATION GRADE)" and "certified" are gone. The overall label is ALL CHECKS
  PASSED, PASSED WITH WARNINGS, AT LEAST ONE CHECK FAILED or NO CHECKS RUN (status NOT_APPLICABLE instead of
  PASS). The methods paragraph states what was checked and found, including warnings and failures; the IRC row
  says that an IRC is the user's statement and is not checked.
- The methods paragraph no longer fills missing metadata with "VASP / ORCA / Gaussian" and "DFT".
- Version strings and citations come from `nebcert.__version__` and the concept DOI (`nebcert/citation.py`).

### Validation (`validation/`)
- `compare_overreact.py`, `eckart_highprecision.py`: Eckart, Wigner and Eyring against overreact and a 40-digit
  mpmath reference on 720 grid points (see README).
- `legacy100_grid.py`: the Eckart factor of 1.0.0 differs from the reference by more than a factor of 2 at 529
  of 720 grid points (median ratio 6.8).
- `make_neb_inputs.py`, `benchmark_neb.py`, `legacy_eval.py`: ORCA 6.1.1 NEB-TS calculations, degraded and failed
  runs, compared with ORCA's full-precision files and with the verdicts of 1.0.0 and 1.1.0.

## 1.1.0 (committed, not released)

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
