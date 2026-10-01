"""
Tests for NEB profile analysis and spline interpolation.
"""

import numpy as np
import pytest
from nebcert.core.neb_profile import calculate_neb_profile_analysis


def test_neb_profile_standard_mep_pass():
    energies = [0.0, 0.2, 0.5, 0.2, -0.3]
    s_coords = [0.0, 0.5, 1.0, 1.5, 2.0]
    forces = [0.0, 0.02, 0.01, 0.02, 0.0]

    res = calculate_neb_profile_analysis(
        energies_ev=energies,
        coordinates_s_ang=s_coords,
        tangent_forces_ev_ang=forces
    )

    assert res.status == "PASS"
    assert res.n_images == 5
    assert np.isclose(res.e_forward_barrier_ev, 0.5, atol=0.01)
    assert np.isclose(res.e_reverse_barrier_ev, 0.8, atol=0.01)
    assert np.isclose(res.delta_e_reaction_ev, -0.3, atol=0.01)
    assert res.is_climbing_image_converged is True


def test_neb_profile_unconverged_forces():
    energies = [0.0, 0.3, 0.6, 0.3, 0.0]
    forces = [0.0, 0.15, 0.25, 0.10, 0.0]  # High forces

    res = calculate_neb_profile_analysis(
        energies_ev=energies,
        tangent_forces_ev_ang=forces,
        force_convergence_threshold=0.05
    )

    assert res.status == "WARNING"
    assert res.is_climbing_image_converged is False


def test_unconverged_band_fails():
    res = calculate_neb_profile_analysis([0.0, 0.2, 0.5, 0.2, -0.3], band_converged=False)
    assert res.status == "FAIL"


def test_unconverged_ts_optimisation_fails():
    res = calculate_neb_profile_analysis([0.0, 0.2, 0.5, 0.2, -0.3], band_converged=True,
                                         ts_energy_rel_ev=0.55, ts_optimisation_converged=False)
    assert res.status == "FAIL"


def test_optimised_ts_energy_replaces_spline_barrier():
    res = calculate_neb_profile_analysis([0.0, 0.2, 0.5, 0.2, -0.3], band_converged=True,
                                         ts_energy_rel_ev=0.55, ts_optimisation_converged=True)
    assert res.status == "PASS"
    assert res.barrier_source == "optimised TS"
    assert np.isclose(res.e_forward_barrier_ev, 0.55)
    assert np.isclose(res.e_reverse_barrier_ev, 0.85)
    assert np.isclose(res.e_forward_barrier_spline_ev, 0.5, atol=0.01)


def test_shallow_dip_is_noise_deep_well_is_intermediate():
    shallow = calculate_neb_profile_analysis([0.0, 0.30, 0.50, 0.49, 0.50, 0.2, -0.3])
    assert shallow.has_intermediate_minimum is False
    deep = calculate_neb_profile_analysis([0.0, 0.30, 0.50, 0.20, 0.45, 0.2, -0.3])
    assert deep.has_intermediate_minimum is True
    assert np.isclose(deep.deepest_intermediate_well_kcal_mol, 0.25 * 23.06054887, rtol=1e-6)
    assert deep.status == "WARNING"


def test_identical_end_points_fail():
    res = calculate_neb_profile_analysis([0.0, 1e-6, 2e-6, 0.0], coordinates_s_ang=[0.0, 0.0002, 0.0004, 0.001],
                                         band_converged=True)
    assert res.status == "FAIL"
    # distances printed as 0.000 for every image (ORCA prints three decimals) must not crash the spline
    res = calculate_neb_profile_analysis([0.0, 0.0, 0.0, 0.0], coordinates_s_ang=[0.0, 0.0, 0.0, 0.001],
                                         band_converged=True)
    assert res.status == "FAIL"
    assert "same structure" in res.diagnostic_message


def test_converged_climbing_image_is_the_barrier():
    e = [0.0, 0.2, 0.5, 0.2, -0.3]
    res = calculate_neb_profile_analysis(e, coordinates_s_ang=[0.0, 0.5, 1.0, 1.5, 2.0], band_converged=True,
                                         climbing_image_index=2)
    assert res.barrier_source == "climbing image"
    assert np.isclose(res.e_forward_barrier_ev, 0.5)
    # without convergence the climbing image is not trusted
    res = calculate_neb_profile_analysis(e, coordinates_s_ang=[0.0, 0.5, 1.0, 1.5, 2.0], climbing_image_index=2)
    assert res.barrier_source == "spline"


def test_band_force_check_skipped_after_converged_ts_optimisation():
    e, f = [0.0, 0.2, 0.5, 0.2, -0.3], [0.0, 0.1, 0.09, 0.1, 0.0]
    assert calculate_neb_profile_analysis(e, tangent_forces_ev_ang=f, band_converged=True).status == "WARNING"
    res = calculate_neb_profile_analysis(e, tangent_forces_ev_ang=f, band_converged=True, ts_energy_rel_ev=0.51,
                                         ts_optimisation_converged=True)
    assert res.status == "PASS"
