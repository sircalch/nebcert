"""
Regression tests for the asymmetric Eckart transmission against an independent reference:
the 1D Schrodinger equation for the Eckart potential solved numerically (transfer matrix,
m = 1 amu), with the imaginary frequency taken from the numerical barrier-top curvature.
"""
import numpy as np
import pytest

from nebcert.core.tunneling import (
    calculate_quantum_tunneling_corrections,
    eckart_transmission_probability,
)

# (V1, V2 [eV], nu [cm^-1], [(E/V1, P_numeric), ...])
REFERENCE = [
    (0.300, 0.300, 1410.0, [(0.50, 2.269e-03), (0.80, 1.144e-01), (1.00, 5.573e-01), (1.20, 9.079e-01)]),
    (0.400, 0.200, 1734.0, [(0.80, 9.410e-02), (1.00, 5.747e-01), (1.20, 9.227e-01)]),
    (0.227, 0.847, 1212.0, [(0.50, 5.460e-03), (0.80, 1.368e-01), (1.00, 5.337e-01), (1.20, 8.773e-01)]),
]


@pytest.mark.parametrize("v1,v2,nu,points", REFERENCE)
def test_transmission_matches_schrodinger_reference(v1, v2, nu, points):
    for frac, p_ref in points:
        p = eckart_transmission_probability(frac * v1, v1, v2, nu)
        # the reference frequencies are rounded to 1 cm^-1, hence the 2% tolerance
        assert p == pytest.approx(p_ref, rel=0.02, abs=5e-5)


def test_transmission_is_zero_below_product_asymptote():
    # Endothermic direction: V1 - V2 = 0.2 eV, no transmission below E = 0.2 eV
    assert eckart_transmission_probability(0.15, 0.4, 0.2, 1500.0) == 0.0


def test_transmission_is_bounded_for_wide_barriers():
    for e in np.linspace(0.01, 3.0, 50):
        p = eckart_transmission_probability(e, 1.5, 1.2, 300.0)
        assert 0.0 <= p <= 1.0


def test_kappa_approaches_wigner_at_high_temperature():
    res = calculate_quantum_tunneling_corrections(1000.0, 0.5, 0.6, temperatures_k=[2000.0])
    pt = res.temperature_points[0]
    assert pt.kappa_eckart == pytest.approx(pt.kappa_wigner, rel=0.05)


def test_kappa_is_moderate_for_typical_h_abstraction():
    # Regression for the former C/delta parameter bug, which gave kappa_Eckart(298 K) ~ 5000
    res = calculate_quantum_tunneling_corrections(1250.0, 0.227, 0.847)
    assert 1.0 < res.kappa_eckart_298 < 20.0
