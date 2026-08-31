"""
Tests for transition state imaginary frequency verification.
"""

import numpy as np
import pytest
from nebcert.core.ts_frequency import verify_ts_frequency_and_irc


def test_ts_frequency_valid_pass():
    freqs = [-850.0, 150.0, 300.0, 600.0, 1200.0, 3000.0]
    res = verify_ts_frequency_and_irc(freqs, irc_confirmed=True)

    assert res.n_imaginary_frequencies == 1
    assert res.is_true_first_order_saddle_point is True
    assert res.is_magnitude_significant is True
    assert res.status == "PASS"


def test_ts_frequency_no_imaginary_fail():
    freqs = [50.0, 150.0, 300.0, 600.0]
    res = verify_ts_frequency_and_irc(freqs)

    assert res.n_imaginary_frequencies == 0
    assert res.is_true_first_order_saddle_point is False
    assert res.status == "FAIL"


def test_ts_frequency_multiple_imaginary_fail():
    freqs = [-500.0, -200.0, 100.0, 300.0]
    res = verify_ts_frequency_and_irc(freqs)

    assert res.n_imaginary_frequencies == 2
    assert res.status == "FAIL"
