"""
Tests for Eyring TST kinetics and Arrhenius fitting.
"""

import numpy as np
import pytest
from nebcert.core.tst_kinetics import calculate_eyring_tst_rates


def test_eyring_tst_kinetics():
    # Barrier = 0.50 eV (~ 11.53 kcal/mol)
    e_act = 0.50
    res = calculate_eyring_tst_rates(e_act)

    assert res.status == "PASS"
    assert res.k_298_s_minus_1 > 0.0
    assert np.isclose(res.e_activation_kcal_mol, 11.53, atol=0.1)
    assert len(res.rate_points) == 9
    assert res.arrhenius_pre_exponential_a_s_minus_1 > 1e11
