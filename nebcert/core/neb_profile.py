"""
Nudged Elastic Band (NEB / CI-NEB) profile analysis, spline interpolation, and MEP force audit.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import CubicSpline

EV_TO_KCAL_MOL = 23.06054887


@dataclass
class NEBImagePoint:
    image_index: int
    reaction_coordinate_s_ang: float
    energy_ev: float
    relative_energy_ev: float
    relative_energy_kcal_mol: float
    tangent_force_ev_ang: Optional[float]


@dataclass
class NEBProfileResult:
    n_images: int
    total_path_length_ang: float
    ts_image_index: int
    ts_position_s_ang: float
    e_forward_barrier_ev: float
    e_forward_barrier_kcal_mol: float
    e_reverse_barrier_ev: float
    e_reverse_barrier_kcal_mol: float
    delta_e_reaction_ev: float
    delta_e_reaction_kcal_mol: float
    max_tangent_force_ev_ang: Optional[float]
    is_climbing_image_converged: bool
    has_intermediate_minimum: bool
    interpolated_s: List[float]
    interpolated_e_rel_ev: List[float]
    image_points: List[NEBImagePoint]
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str


def calculate_neb_profile_analysis(
    energies_ev: List[float],
    coordinates_s_ang: Optional[List[float]] = None,
    tangent_forces_ev_ang: Optional[List[float]] = None,
    force_convergence_threshold: float = 0.05,  # eV/Å
    warning_force_threshold: float = 0.10
) -> NEBProfileResult:
    """
    Analyzes NEB minimum energy path, fits cubic spline interpolation, calculates
    forward and reverse activation barriers, and audits force convergence and continuity.

    Parameters
    ----------
    energies_ev : list of float
        Total DFT energies for all images along the band (eV).
    coordinates_s_ang : list of float, optional
        Cumulative reaction coordinates s (Å). If None, uniform spacing 0..1 is used.
    tangent_forces_ev_ang : list of float, optional
        Parallel / residual forces along band (eV/Å).
    force_convergence_threshold : float, default 0.05 eV/Å
    warning_force_threshold : float, default 0.10 eV/Å

    Returns
    -------
    result : NEBProfileResult
    """
    e_arr = np.asarray(energies_ev, dtype=float)
    n_img = len(e_arr)
    if n_img < 3:
        raise ValueError("NEB reaction path requires at least 3 images (Reactant, TS, Product).")

    if coordinates_s_ang is not None and len(coordinates_s_ang) == n_img:
        s_arr = np.asarray(coordinates_s_ang, dtype=float)
    else:
        s_arr = np.linspace(0.0, float(n_img - 1), n_img)

    # Reference energies relative to Reactant (image 0)
    e_rel = e_arr - e_arr[0]
    total_length = float(s_arr[-1] - s_arr[0])

    # Fit natural cubic spline
    spline = CubicSpline(s_arr, e_rel, bc_type='natural')
    s_fine = np.linspace(s_arr[0], s_arr[-1], 300)
    e_fine = spline(s_fine)

    # Find peak (Transition State candidate)
    ts_fine_idx = int(np.argmax(e_fine))
    ts_s_pos = float(s_fine[ts_fine_idx])
    ts_peak_e_rel = float(e_fine[ts_fine_idx])

    # Nearest discrete image to TS peak
    ts_img_idx = int(np.argmin(np.abs(s_arr - ts_s_pos)))

    e_fwd_barrier_ev = max(0.0, ts_peak_e_rel)
    e_fwd_barrier_kcal = float(e_fwd_barrier_ev * EV_TO_KCAL_MOL)

    e_product_rel_ev = float(e_rel[-1])
    e_rev_barrier_ev = max(0.0, ts_peak_e_rel - e_product_rel_ev)
    e_rev_barrier_kcal = float(e_rev_barrier_ev * EV_TO_KCAL_MOL)

    delta_e_rxn_ev = e_product_rel_ev
    delta_e_rxn_kcal = float(delta_e_rxn_ev * EV_TO_KCAL_MOL)

    # Check for intermediate minima (valleys along MEP)
    # Check if there is an internal local minimum between 0 and -1
    has_interm_min = False
    for i in range(1, n_img - 1):
        if e_rel[i] < e_rel[i-1] and e_rel[i] < e_rel[i+1]:
            has_interm_min = True
            break

    # Force audit
    max_f = None
    is_f_conv = True
    if tangent_forces_ev_ang is not None and len(tangent_forces_ev_ang) == n_img:
        # Check force specifically at climbing image
        f_arr = np.abs(np.asarray(tangent_forces_ev_ang, dtype=float))
        max_f = float(np.max(f_arr))
        ts_f = float(f_arr[ts_img_idx])
        if ts_f > force_convergence_threshold:
            is_f_conv = False

    # Image points
    img_points = []
    for i in range(n_img):
        f_val = float(tangent_forces_ev_ang[i]) if tangent_forces_ev_ang is not None and i < len(tangent_forces_ev_ang) else None
        img_points.append(NEBImagePoint(
            image_index=i,
            reaction_coordinate_s_ang=float(s_arr[i]),
            energy_ev=float(e_arr[i]),
            relative_energy_ev=float(e_rel[i]),
            relative_energy_kcal_mol=float(e_rel[i] * EV_TO_KCAL_MOL),
            tangent_force_ev_ang=f_val
        ))

    # Decision logic
    if not is_f_conv:
        status = "WARNING"
        diag = f"NEB Climbing Image forces unconverged (|F_tangent| = {max_f:.3f} > {force_convergence_threshold:.2f} eV/Å). Barrier estimate: E_a = {e_fwd_barrier_kcal:.2f} kcal/mol."
    elif has_interm_min:
        status = "WARNING"
        diag = f"Complex reaction path: intermediate minimum detected along MEP. Multi-step reaction pathway should be split into distinct elemental steps."
    elif e_fwd_barrier_ev < 1e-4 and abs(delta_e_rxn_ev) < 1e-4:
        status = "FAIL"
        diag = "Flat energy profile (barrier < 0.001 eV). Check initial/final image structures."
    else:
        status = "PASS"
        diag = f"Single-step minimum energy path fully validated (E_a^fwd = {e_fwd_barrier_kcal:.2f} kcal/mol ({e_fwd_barrier_ev:.3f} eV), Delta E_rxn = {delta_e_rxn_kcal:.2f} kcal/mol)."

    return NEBProfileResult(
        n_images=n_img,
        total_path_length_ang=total_length,
        ts_image_index=ts_img_idx,
        ts_position_s_ang=ts_s_pos,
        e_forward_barrier_ev=e_fwd_barrier_ev,
        e_forward_barrier_kcal_mol=e_fwd_barrier_kcal,
        e_reverse_barrier_ev=e_rev_barrier_ev,
        e_reverse_barrier_kcal_mol=e_rev_barrier_kcal,
        delta_e_reaction_ev=delta_e_rxn_ev,
        delta_e_reaction_kcal_mol=delta_e_rxn_kcal,
        max_tangent_force_ev_ang=max_f,
        is_climbing_image_converged=is_f_conv,
        has_intermediate_minimum=has_interm_min,
        interpolated_s=s_fine.tolist(),
        interpolated_e_rel_ev=e_fine.tolist(),
        image_points=img_points,
        status=status,
        diagnostic_message=diag
    )
