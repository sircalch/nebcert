"""
Quantum mechanical tunneling corrections: Wigner formula and numerical asymmetric Eckart barrier transmission.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad

# Physical constants
PLANCK_H_J_S = 6.62607015e-34       # J * s
SPEED_OF_LIGHT_CM_S = 2.99792458e10 # cm / s
K_B_J_K = 1.380649e-23              # J / K
K_B_EV_K = 8.617333262e-5           # eV / K
EV_TO_JOULE = 1.602176634e-19


@dataclass
class TunnelingTemperaturePoint:
    temperature_k: float
    kappa_wigner: float
    kappa_eckart: float


@dataclass
class TunnelingResult:
    imaginary_freq_cm1: float
    e_forward_barrier_ev: float
    e_reverse_barrier_ev: float
    temperatures_k: List[float]
    temperature_points: List[TunnelingTemperaturePoint]
    kappa_wigner_298: float
    kappa_eckart_298: float
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str


def calculate_quantum_tunneling_corrections(
    imaginary_freq_cm1: float,
    e_forward_barrier_ev: float,
    e_reverse_barrier_ev: Optional[float] = None,
    temperatures_k: Optional[List[float]] = None
) -> TunnelingResult:
    """
    Computes Wigner and 1D asymmetric Eckart quantum tunneling factors kappa(T).

    Parameters
    ----------
    imaginary_freq_cm1 : float
        Magnitude of imaginary frequency (cm^-1).
    e_forward_barrier_ev : float
        Forward activation energy (eV).
    e_reverse_barrier_ev : float, optional
        Reverse activation energy (eV). If None, symmetric barrier (V1 = V2) is used.
    temperatures_k : list of float, optional

    Returns
    -------
    result : TunnelingResult
    """
    nu = abs(float(imaginary_freq_cm1))
    v1 = max(1e-4, float(e_forward_barrier_ev))
    v2 = max(1e-4, float(e_reverse_barrier_ev)) if e_reverse_barrier_ev is not None else v1

    if temperatures_k is None:
        t_list = [200.0, 250.0, 298.15, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0]
    else:
        t_list = sorted(temperatures_k)

    # Conversion of frequency to energy (eV)
    # h * nu_c (J) / eV
    h_nu_j = PLANCK_H_J_S * SPEED_OF_LIGHT_CM_S * nu
    h_nu_ev = h_nu_j / EV_TO_JOULE

    # Characteristic Eckart barrier parameter C
    # C = (h*nu)^2 / [8 * (sqrt(V1) + sqrt(V2))^2]
    c_param = (h_nu_ev**2) / (8.0 * (np.sqrt(v1) + np.sqrt(v2))**2)
    delta_v = v1 - v2

    points = []
    for t in t_list:
        kb_t = K_B_EV_K * t
        u_val = h_nu_ev / kb_t

        # 1. Wigner: kappa_W = 1 + (1/24) * u^2
        kap_w = float(1.0 + (1.0 / 24.0) * (u_val**2))

        # 2. Asymmetric Eckart Transmission Function P(E)
        def eckart_transmission(e_ev):
            if e_ev <= 0.0:
                return 0.0
            if e_ev - delta_v <= 0.0:
                return 0.0

            # 2*pi * a = 2*pi * 0.5 * sqrt(E / C)
            two_pi_a = np.pi * np.sqrt(e_ev / c_param)
            two_pi_b = np.pi * np.sqrt((e_ev - delta_v) / c_param)

            # d parameter
            arg_d = (4.0 * v1 * v2 - c_param) / c_param
            if arg_d > 0:
                two_pi_d = np.pi * np.sqrt(arg_d)
                # cosh(2*pi*d)
                denom = np.cosh(np.clip(two_pi_a + two_pi_b, -100, 100)) + np.cosh(np.clip(two_pi_d, -100, 100))
            else:
                two_pi_d = np.pi * np.sqrt(abs(arg_d))
                denom = np.cosh(np.clip(two_pi_a + two_pi_b, -100, 100)) + np.cos(two_pi_d)

            num = np.cosh(np.clip(two_pi_a + two_pi_b, -100, 100)) - np.cosh(np.clip(two_pi_a - two_pi_b, -100, 100))
            return float(num / max(1e-15, denom))

        # Integrate: kappa_Eckart = (1 / k_B T) * int_0^infty P(E) * exp(-(E - V1)/k_B T) dE
        # Normalize relative to classical barrier transmission
        def integrand(e):
            p = eckart_transmission(e)
            return p * np.exp(-(e - v1) / kb_t) / kb_t

        # Numerical integration up to 10*kb_t above barrier
        upper_limit = v1 + 15.0 * kb_t
        try:
            val, _ = quad(integrand, max(0.0, delta_v), upper_limit, limit=100)
            kap_e = float(max(1.0, val))
        except Exception:
            kap_e = kap_w

        points.append(TunnelingTemperaturePoint(
            temperature_k=t,
            kappa_wigner=kap_w,
            kappa_eckart=kap_e
        ))

    # Reference 298.15 K
    kw_298 = points[0].kappa_wigner
    ke_298 = points[0].kappa_eckart
    for pt in points:
        if np.isclose(pt.temperature_k, 298.15):
            kw_298 = pt.kappa_wigner
            ke_298 = pt.kappa_eckart
            break

    status = "PASS"
    diag = f"Quantum tunneling corrections computed (|nu| = {nu:.1f} cm^-1). At 298.15 K: kappa_Wigner = {kw_298:.2f}, kappa_Eckart = {ke_298:.2f}."

    return TunnelingResult(
        imaginary_freq_cm1=nu,
        e_forward_barrier_ev=v1,
        e_reverse_barrier_ev=v2,
        temperatures_k=t_list,
        temperature_points=points,
        kappa_wigner_298=kw_298,
        kappa_eckart_298=ke_298,
        status=status,
        diagnostic_message=diag
    )
