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
    Counts the imaginary frequencies of a candidate transition state. PASS requires exactly one with
    |nu| >= min_significant_freq_cm1 and no smaller imaginary modes.

    Parameters
    ----------
    frequencies_cm1 : list of float
        Vibrational harmonic frequencies in cm^-1 (negative values represent imaginary frequencies).
    min_significant_freq_cm1 : float, default 50.0 cm^-1
    irc_confirmed : bool, optional
        Whether an IRC calculation connected the intended reactant and product minima, as reported by the
        user; NEBCert does not read or check the IRC itself.

    Returns
    -------
    result : TSFrequencyResult
    """
    freqs = np.asarray(frequencies_cm1, dtype=float)
    # Imaginary modes are printed as negative numbers. Modes between -min_significant and
    # -noise_floor are counted separately: they are usually numerical noise (loose optimisation,
    # integration grid) rather than a second reaction coordinate, but cannot be ignored silently.
    noise_floor = 1.0
    sig_mask = freqs <= -min_significant_freq_cm1
    small_mask = (freqs < -noise_floor) & ~sig_mask
    n_sig = int(np.sum(sig_mask))
    n_small = int(np.sum(small_mask))
    n_imag = n_sig + n_small

    imag_val = None
    if n_imag >= 1:
        imag_val = float(np.min(freqs[freqs < -noise_floor]))
    is_sig = n_sig >= 1
    is_first_order = (n_sig == 1 and n_small == 0)
    small_str = ", ".join(f"{f:.1f}" for f in sorted(freqs[small_mask]))

    if n_imag == 0:
        status = "FAIL"
        diag = "No imaginary frequency: the structure is a minimum, not a transition state."
    elif n_sig > 1:
        status = "FAIL"
        diag = (f"{n_sig} imaginary frequencies of magnitude >= {min_significant_freq_cm1:.0f} cm^-1 "
                f"(higher-order saddle point).")
    elif n_sig == 0:
        status = "WARNING"
        diag = (f"Only small imaginary frequencies ({small_str} cm^-1, |nu| < {min_significant_freq_cm1:.0f}): "
                f"no clear reaction-coordinate mode.")
    elif n_small > 0:
        status = "WARNING"
        diag = (f"One imaginary frequency of {imag_val:.1f} cm^-1 plus {n_small} small one(s) ({small_str} cm^-1); "
                f"tighten the optimisation or the grid and recompute the frequencies.")
    elif irc_confirmed is False:
        status = "WARNING"
        diag = f"One imaginary frequency ({imag_val:.1f} cm^-1), but the IRC did not connect the intended minima."
    else:
        status = "PASS"
        irc_str = (" IRC connectivity reported by the user." if irc_confirmed is True else
                   " Whether this mode connects the intended minima is not checked (no IRC).")
        diag = f"One imaginary frequency ({imag_val:.1f} cm^-1).{irc_str}"

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
