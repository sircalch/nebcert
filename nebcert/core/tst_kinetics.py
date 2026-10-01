"""
Transition State Theory (TST / Eyring-Polanyi) rate constants and Arrhenius parameter fitting.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np

# Physical constants
K_B_EV_K = 8.617333262e-5      # eV / K
K_B_J_K = 1.380649e-23         # J / K
PLANCK_H_J_S = 6.62607015e-34  # J * s
GAS_CONST_KCAL = 1.98720425864083e-3  # kcal / (mol * K)
EV_TO_KCAL_MOL = 23.06054887


@dataclass
class TSTTemperaturePoint:
    temperature_k: float
    inv_temperature_1_k: float
    k_tst_s_minus_1: float
    log_k_tst: float
    half_life_seconds: Optional[float]


@dataclass
class TSTKineticsResult:
    e_activation_ev: float
    e_activation_kcal_mol: float
    temperatures_k: List[float]
    rate_points: List[TSTTemperaturePoint]
    k_298_s_minus_1: float
    half_life_298_s: float
    arrhenius_pre_exponential_a_s_minus_1: float
    arrhenius_e_activation_kcal_mol: float
    status: str  # 'PASS' or 'NOT_APPLICABLE'
    diagnostic_message: str


def calculate_eyring_tst_rates(
    e_activation_ev: float,
    temperatures_k: Optional[List[float]] = None,
    barrier_type: str = "electronic"
) -> TSTKineticsResult:
    """
    Computes Eyring-Polanyi transition state theory rate constants k(T) across temperatures
    and extracts Arrhenius parameters A and E_a.

    Parameters
    ----------
    e_activation_ev : float
        Barrier inserted in k = (k_B T / h) exp(-E / k_B T), in eV.
    temperatures_k : list of float, optional
        List of temperatures in Kelvin (default: [200, 250, 298.15, 350, 400, 500, 600, 800, 1000]).
    barrier_type : {'electronic', 'gibbs'}
        'gibbs' when the barrier is a Gibbs free energy of activation (Delta G double-dagger): the
        result is the Eyring rate constant (PASS). 'electronic' (the default, and what an NEB band
        gives) means an electronic energy difference without zero-point, thermal or entropic terms:
        the number is then not a TST rate constant and the status is NOT_APPLICABLE (reported, but not
        counted in the overall status).

    Returns
    -------
    result : TSTKineticsResult
    """
    if temperatures_k is None:
        t_list = [200.0, 250.0, 298.15, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0]
    else:
        t_list = sorted(temperatures_k)

    e_act_kcal = float(e_activation_ev * EV_TO_KCAL_MOL)
    rate_pts = []

    for t in t_list:
        # Pre-exponential frequency factor k_B * T / h
        freq_factor = (K_B_J_K * t) / PLANCK_H_J_S  # ~ 6.21e12 s^-1 at 298.15 K
        
        # Boltzmann exponent
        exp_term = np.exp(-e_activation_ev / (K_B_EV_K * t))
        k_val = float(freq_factor * exp_term)
        log_k = float(np.log(max(1e-300, k_val)))
        
        # First order half life t_1/2 = ln(2) / k
        t_half = float(np.log(2.0) / k_val) if k_val > 1e-100 else 1e30

        rate_pts.append(TSTTemperaturePoint(
            temperature_k=t,
            inv_temperature_1_k=float(1000.0 / t),  # 1000 / T for Arrhenius plot
            k_tst_s_minus_1=k_val,
            log_k_tst=log_k,
            half_life_seconds=t_half
        ))

    # Arrhenius linear fit: ln(k) = ln(A) - (E_a / R) * (1 / T)
    inv_t = np.array([1.0 / pt.temperature_k for pt in rate_pts])
    ln_k = np.array([pt.log_k_tst for pt in rate_pts])
    poly = np.polyfit(inv_t, ln_k, 1)
    slope = poly[0]
    intercept = poly[1]

    arr_a = float(np.exp(min(700.0, intercept)))
    arr_ea_kcal = float(-slope * GAS_CONST_KCAL)

    # 298.15 K reference values
    k_298 = None
    half_298 = None
    for pt in rate_pts:
        if np.isclose(pt.temperature_k, 298.15):
            k_298 = pt.k_tst_s_minus_1
            half_298 = pt.half_life_seconds
            break
    if k_298 is None:
        k_298 = rate_pts[0].k_tst_s_minus_1
        half_298 = rate_pts[0].half_life_seconds

    if barrier_type not in ("electronic", "gibbs"):
        raise ValueError("barrier_type must be 'electronic' or 'gibbs'")
    summary = (f"k_298 = {k_298:.2e} s^-1 over {len(t_list)} temperatures (Arrhenius fit A = {arr_a:.2e} s^-1, "
               f"E_a = {arr_ea_kcal:.2f} kcal/mol)")
    if barrier_type == "gibbs":
        status = "PASS"
        diag = f"Eyring rate constant from a Gibbs barrier of {e_act_kcal:.2f} kcal/mol: {summary}."
    else:
        status = "NOT_APPLICABLE"
        diag = (f"Eyring expression evaluated with an electronic barrier of {e_act_kcal:.2f} kcal/mol "
                f"(no zero-point, thermal or entropic terms): {summary}. This is not a TST rate constant; "
                f"supply a Gibbs barrier for one.")

    return TSTKineticsResult(
        e_activation_ev=e_activation_ev,
        e_activation_kcal_mol=e_act_kcal,
        temperatures_k=t_list,
        rate_points=rate_pts,
        k_298_s_minus_1=k_298,
        half_life_298_s=half_298,
        arrhenius_pre_exponential_a_s_minus_1=arr_a,
        arrhenius_e_activation_kcal_mol=arr_ea_kcal,
        status=status,
        diagnostic_message=diag
    )
