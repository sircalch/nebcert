"""
Publication-ready vector figures for NEB reaction profiles, Arrhenius kinetics, and tunneling factors.
"""

from typing import List, Optional
import os
import numpy as np
import matplotlib.pyplot as plt
from nebcert.core.scoring import ReactionPathwayReport

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'lines.linewidth': 2.0,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})

EV_TO_KCAL_MOL = 23.06054887


def generate_nebcert_figures(
    report: ReactionPathwayReport,
    output_dir: str,
    formats: List[str] = ("png", "svg", "pdf")
) -> List[str]:
    """
    Generates publication figures: NEB minimum energy path, Arrhenius plot, and tunneling factors.

    Parameters
    ----------
    report : ReactionPathwayReport
    output_dir : str
    formats : list of str

    Returns
    -------
    saved_paths : list of str
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = []

    # 1. NEB Minimum Energy Path Profile
    if report.neb_profile is not None:
        neb = report.neb_profile
        s_pts = [pt.reaction_coordinate_s_ang for pt in neb.image_points]
        e_pts = [pt.relative_energy_kcal_mol for pt in neb.image_points]
        
        s_fine = np.asarray(neb.interpolated_s)
        e_fine_kcal = np.asarray(neb.interpolated_e_rel_ev) * EV_TO_KCAL_MOL

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(s_fine, e_fine_kcal, color="#0284c7", linewidth=2.2, label="Spline Interpolated MEP")
        ax.scatter(s_pts, e_pts, color="#0f172a", s=60, zorder=5, label="NEB Images")

        # Highlight TS
        ts_s = neb.ts_position_s_ang
        ts_e = neb.e_forward_barrier_kcal_mol
        ax.scatter([ts_s], [ts_e], color="#dc2626", s=120, marker="*", zorder=6, label=rf"TS ($E_a^{{\mathrm{{fwd}}}} = {ts_e:.2f}\ \mathrm{{kcal/mol}}$)")

        # Barrier annotation
        ax.annotate(
            rf"$E_a^{{\mathrm{{fwd}}}} = {neb.e_forward_barrier_kcal_mol:.2f}\ \mathrm{{kcal/mol}}$",
            xy=(ts_s, ts_e),
            xytext=(ts_s, ts_e * 1.08),
            ha="center",
            fontweight="bold",
            color="#dc2626",
            fontsize=10
        )

        ax.set_xlabel(r"Reaction Coordinate $s$ ($\mathrm{\AA}$)")
        ax.set_ylabel(r"Relative Energy $\Delta E$ (kcal/mol)")
        rxn_name = report.metadata.get("reaction", "Reaction Pathway")
        ax.set_title(f"NEB Minimum Energy Path — {rxn_name}")
        ax.grid(True)
        ax.legend(loc="upper right", frameon=True, fontsize=9)

        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"nebcert_neb_mep_profile.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 2. Arrhenius Kinetics Plot ln(k) vs 1000/T
    if report.tst_kinetics is not None:
        tst = report.tst_kinetics
        inv_t = [pt.inv_temperature_1_k for pt in tst.rate_points]
        ln_k = [pt.log_k_tst for pt in tst.rate_points]

        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.plot(inv_t, ln_k, marker="o", markersize=6, color="#16a34a", linewidth=2.0, label=r"Classical TST $k_{\mathrm{TST}}(T)$")

        # Overlay tunneling-corrected rate if available
        if report.tunneling is not None:
            tun = report.tunneling
            t_pts = tun.temperature_points
            ln_k_tun = []
            for i, pt in enumerate(tst.rate_points):
                if i < len(t_pts):
                    kap = t_pts[i].kappa_eckart
                    ln_k_tun.append(pt.log_k_tst + np.log(kap))
                else:
                    ln_k_tun.append(pt.log_k_tst)
            ax.plot(inv_t, ln_k_tun, marker="s", markersize=6, color="#dc2626", linestyle="--", linewidth=2.0, label=r"Tunneling Corrected $k_{\mathrm{Eckart}}(T)$")

        ax.set_xlabel(r"$1000 / T\ (\mathrm{K}^{-1})$")
        ax.set_ylabel(r"$\ln(k\ /\ \mathrm{s}^{-1})$")
        ax.set_title(rf"Arrhenius Kinetics — $A = {tst.arrhenius_pre_exponential_a_s_minus_1:.2e}\ \mathrm{{s}}^{{-1}}$, $E_a = {tst.arrhenius_e_activation_kcal_mol:.2f}\ \mathrm{{kcal/mol}}$")
        ax.grid(True)
        ax.legend(loc="upper right", frameon=True, fontsize=9)

        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"nebcert_arrhenius_kinetics.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 3. Quantum Tunneling Factor Curve kappa(T)
    if report.tunneling is not None:
        tun = report.tunneling
        t_vals = [pt.temperature_k for pt in tun.temperature_points]
        kw_vals = [pt.kappa_wigner for pt in tun.temperature_points]
        ke_vals = [pt.kappa_eckart for pt in tun.temperature_points]

        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ax.plot(t_vals, ke_vals, color="#0284c7", linewidth=2.2, label=r"Eckart Transmission $\kappa_{\mathrm{Eckart}}(T)$")
        ax.plot(t_vals, kw_vals, color="#f59e0b", linestyle="--", linewidth=2.0, label=r"Wigner Factor $\kappa_{\mathrm{Wigner}}(T)$")
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=1.0, label="Classical Limit (No Tunneling)")

        ax.set_xlabel("Temperature $T$ (K)")
        ax.set_ylabel(r"Tunneling Transmission Factor $\kappa(T)$")
        ax.set_title(rf"Quantum Tunneling Correction — $|\nu^\ddagger| = {tun.imaginary_freq_cm1:.1f}\ \mathrm{{cm}}^{{-1}}$")
        ax.grid(True)
        ax.legend(loc="upper right", frameon=True, fontsize=9)

        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"nebcert_tunneling_factors.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    return saved_files
