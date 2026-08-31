"""
Tests for Wigner and Eckart quantum tunneling corrections.
"""

import numpy as np
import pytest
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections


def test_quantum_tunneling_corrections():
    res = calculate_quantum_tunneling_corrections(
        imaginary_freq_cm1=1000.0,
        e_forward_barrier_ev=0.40,
        e_reverse_barrier_ev=0.60
    )

    assert res.status == "PASS"
    assert res.kappa_wigner_298 > 1.0
    assert res.kappa_eckart_298 >= 1.0
    assert len(res.temperature_points) == 9
