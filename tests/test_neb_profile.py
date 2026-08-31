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
