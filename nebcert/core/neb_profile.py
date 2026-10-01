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
    e_forward_barrier_spline_ev: Optional[float] = None
    barrier_source: str = "spline"
    band_converged: Optional[bool] = None
    ts_optimisation_converged: Optional[bool] = None
    deepest_intermediate_well_kcal_mol: float = 0.0


def calculate_neb_profile_analysis(
    energies_ev: List[float],
    coordinates_s_ang: Optional[List[float]] = None,
    tangent_forces_ev_ang: Optional[List[float]] = None,
    force_convergence_threshold: float = 0.05,  # eV/Å
    warning_force_threshold: float = 0.10,
    band_converged: Optional[bool] = None,
    ts_energy_rel_ev: Optional[float] = None,
    ts_optimisation_converged: Optional[bool] = None,
    min_well_depth_kcal_mol: float = 0.5,
    min_path_length_ang: float = 0.05,
    climbing_image_index: Optional[int] = None
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
        Limit on the force at the image nearest the barrier.
    warning_force_threshold : float, default 0.10 eV/Å
        Unused; kept for backward compatibility.
    band_converged : bool, optional
        Convergence flag of the band as printed by the program (ORCA: "THE NEB OPTIMIZATION HAS
        CONVERGED"). False gives FAIL: an unconverged band does not locate the barrier.
    ts_energy_rel_ev : float, optional
        Energy of an optimised transition state (e.g. ORCA NEB-TS) relative to image 0. When given it
        replaces the spline maximum as the forward barrier; the spline value is kept for comparison.
    ts_optimisation_converged : bool, optional
        False gives FAIL (the TS energy is then not a saddle-point energy). When True, the force check on the
        band is skipped: NEB-TS converges the band only loosely by design and the TS optimisation refines it.
    climbing_image_index : int, optional
        Index of the climbing image. When the band converged and no optimised TS is given, the climbing-image
        energy is the barrier (a converged climbing image sits at the saddle point; a spline through a coarse
        band does not).
    min_well_depth_kcal_mol : float, default 0.5
        An interior image lower than both neighbours by more than this is reported as an intermediate
        minimum; shallower dips are treated as noise of the band.
    min_path_length_ang : float, default 0.05
        When coordinates_s_ang are distances along the path (Å), a shorter path means that both end points
        are the same structure (e.g. one end point slid into the other during pre-optimisation): FAIL.

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
    # Distances printed with few decimals can repeat (e.g. a band whose end points coincide); the spline
    # then needs a strictly increasing abscissa, so it falls back to the image index.
    if not np.all(np.diff(s_arr) > 0):
        s_arr = np.linspace(0.0, float(n_img - 1), n_img)

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

    spline_barrier_ev = max(0.0, ts_peak_e_rel)
    barrier_source = "spline"
    if ts_energy_rel_ev is not None:
        ts_peak_e_rel = float(ts_energy_rel_ev)
        barrier_source = "optimised TS"
    elif climbing_image_index is not None and band_converged and 0 < climbing_image_index < n_img - 1:
        ts_peak_e_rel = float(e_rel[climbing_image_index])
        barrier_source = "climbing image"
    e_fwd_barrier_ev = max(0.0, ts_peak_e_rel)
    e_fwd_barrier_kcal = float(e_fwd_barrier_ev * EV_TO_KCAL_MOL)

    e_product_rel_ev = float(e_rel[-1])
    e_rev_barrier_ev = max(0.0, ts_peak_e_rel - e_product_rel_ev)
    e_rev_barrier_kcal = float(e_rev_barrier_ev * EV_TO_KCAL_MOL)

    delta_e_rxn_ev = e_product_rel_ev
    delta_e_rxn_kcal = float(delta_e_rxn_ev * EV_TO_KCAL_MOL)

    # Check for intermediate minima (valleys along MEP)
    # Check if there is an internal local minimum between 0 and -1
    well = 0.0
    for i in range(1, n_img - 1):
        if e_rel[i] < e_rel[i-1] and e_rel[i] < e_rel[i+1]:
            # depth below the lower of the two barriers that bound it
            depth = min(np.max(e_rel[:i]), np.max(e_rel[i+1:])) - e_rel[i]
            well = max(well, float(depth * EV_TO_KCAL_MOL))
    has_interm_min = well > min_well_depth_kcal_mol

    # Force audit
    max_f = None
    ts_f = None
    is_f_conv = True
    if tangent_forces_ev_ang is not None and len(tangent_forces_ev_ang) == n_img:
        # Check force specifically at climbing image
        f_arr = np.abs(np.asarray(tangent_forces_ev_ang, dtype=float))
        max_f = float(np.max(f_arr))
        ts_f = float(f_arr[ts_img_idx])
        if ts_f > force_convergence_threshold and ts_optimisation_converged is not True:
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

    # Decision logic (first matching condition wins)
    bar = f"E_a^fwd = {e_fwd_barrier_kcal:.2f} kcal/mol ({barrier_source})"
    if coordinates_s_ang is not None and total_length < min_path_length_ang:
        status = "FAIL"
        diag = (f"The path is {total_length:.4f} Å long: both end points are the same structure, so the band "
                f"describes no reaction. Check the end-point optimisations.")
    elif band_converged is False:
        status = "FAIL"
        diag = f"The band did not converge; the highest image does not locate the barrier ({bar})."
    elif ts_optimisation_converged is False:
        status = "FAIL"
        diag = f"The transition-state optimisation did not converge ({bar})."
    elif e_fwd_barrier_ev < 1e-4 and abs(delta_e_rxn_ev) < 1e-4:
        status = "FAIL"
        diag = "Flat energy profile (barrier and reaction energy < 1e-4 eV). Check the end points."
    elif not is_f_conv:
        status = "WARNING"
        diag = (f"Force at the image nearest the barrier is {ts_f:.3f} eV/Å > {force_convergence_threshold:.2f} "
                f"eV/Å ({bar}).")
    elif has_interm_min:
        status = "WARNING"
        diag = (f"Intermediate minimum {well:.2f} kcal/mol deep along the band: the path has more than one "
                f"step, which should be located separately ({bar}).")
    else:
        status = "PASS"
        conv = " Band convergence not reported." if band_converged is None else ""
        diag = (f"No problem found in the band: {bar}, Delta E_rxn = {delta_e_rxn_kcal:.2f} kcal/mol, "
                f"no intermediate minimum deeper than {min_well_depth_kcal_mol:.1f} kcal/mol.{conv}")

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
        diagnostic_message=diag,
        e_forward_barrier_spline_ev=float(spline_barrier_ev),
        barrier_source=barrier_source,
        band_converged=band_converged,
        ts_optimisation_converged=ts_optimisation_converged,
        deepest_intermediate_well_kcal_mol=well
    )
