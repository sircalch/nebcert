"""
Transition State (TS) imaginary frequency verification and IRC connectivity audit.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class TSFrequencyResult:
    n_imaginary_frequencies: int
    imaginary_frequency_cm1: Optional[float]
    all_frequencies_cm1: List[float]
    is_true_first_order_saddle_point: bool
    is_magnitude_significant: bool
    irc_confirmed: Optional[bool]
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str


def verify_ts_frequency_and_irc(
    frequencies_cm1: List[float],
    min_significant_freq_cm1: float = 50.0,
    irc_confirmed: Optional[bool] = None
) -> TSFrequencyResult:
    """
    Verifies that a candidate transition state has strictly 1 imaginary frequency
    with significant magnitude along the reaction coordinate, and audits IRC validation.

    Parameters
    ----------
    frequencies_cm1 : list of float
        Vibrational harmonic frequencies in cm^-1 (negative values represent imaginary frequencies).
    min_significant_freq_cm1 : float, default 50.0 cm^-1
    irc_confirmed : bool, optional
        Whether an IRC calculation successfully verified connection to reactant and product minima.

    Returns
    -------
    result : TSFrequencyResult
    """
    freqs = np.asarray(frequencies_cm1, dtype=float)
    imag_mask = freqs < -1.0  # Threshold to exclude numerical zero modes
    n_imag = int(np.sum(imag_mask))

    imag_val = None
    is_sig = False
    if n_imag >= 1:
        imag_val = float(np.min(freqs[imag_mask]))
        if abs(imag_val) >= min_significant_freq_cm1:
            is_sig = True

    is_first_order = (n_imag == 1)

    # Certification decision
    if n_imag == 0:
        status = "FAIL"
        diag = "No imaginary frequencies detected (0 imaginary modes). Structure is a local minimum, NOT a transition state."
    elif n_imag > 1:
        status = "FAIL"
        diag = f"Higher-order saddle point detected ({n_imag} imaginary frequencies). Structure is an artifact and not a true transition state."
    elif not is_sig:
        status = "WARNING"
        diag = f"Weak imaginary frequency (|nu_imag| = {abs(imag_val):.1f} < {min_significant_freq_cm1:.1f} cm^-1). Likely an unphysical torsional or soft lattice mode rather than bond breaking/formation."
    elif irc_confirmed is False:
        status = "WARNING"
        diag = f"Single imaginary frequency verified ({imag_val:.1f} cm^-1), but IRC path failed to connect intended reactants and products."
    else:
        status = "PASS"
        irc_str = " (IRC connectivity confirmed)" if irc_confirmed is True else ""
        diag = f"Strict 1st-order saddle point certified with 1 significant imaginary frequency (nu_imag = {imag_val:.1f} cm^-1){irc_str}."

    return TSFrequencyResult(
        n_imaginary_frequencies=n_imag,
        imaginary_frequency_cm1=imag_val,
        all_frequencies_cm1=freqs.tolist(),
        is_true_first_order_saddle_point=is_first_order,
        is_magnitude_significant=is_sig,
        irc_confirmed=irc_confirmed,
        status=status,
        diagnostic_message=diag
    )
