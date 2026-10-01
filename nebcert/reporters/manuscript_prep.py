"""
Manuscript Methods snippet, summary tables (CSV, LaTeX), and BibTeX citations for NEBCert.
"""

from typing import Dict, Any, Optional
import os
import pandas as pd
from nebcert import __version__
from nebcert.citation import BIBTEX
from nebcert.core.scoring import ReactionPathwayReport


def generate_nebcert_manuscript_assets(
    report: ReactionPathwayReport,
    output_dir: str
) -> Dict[str, str]:
    """
    Generates manuscript Methods paragraph, summary CSV/LaTeX tables, and BibTeX citations.

    Parameters
    ----------
    report : ReactionPathwayReport
    output_dir : str

    Returns
    -------
    paths : dict
    """
    os.makedirs(output_dir, exist_ok=True)
    generated = {}

    rows = []
    meta = report.metadata

    rows.append({"Parameter": "Reaction", "Value": f"{meta.get('reaction', 'Chemical Reaction')} ({meta.get('software', 'not given')})", "Status": ""})
    rows.append({"Parameter": "Level of theory", "Value": f"{meta.get('functional', 'not given')}", "Status": ""})

    if report.neb_profile:
        neb = report.neb_profile
        rows.append({"Parameter": "Forward Activation Barrier (E_a^fwd)", "Value": f"{neb.e_forward_barrier_kcal_mol:.2f} kcal/mol ({neb.e_forward_barrier_ev:.3f} eV, {neb.barrier_source})", "Status": neb.status})
        rows.append({"Parameter": "Reverse Activation Barrier (E_a^rev)", "Value": f"{neb.e_reverse_barrier_kcal_mol:.2f} kcal/mol ({neb.e_reverse_barrier_ev:.3f} eV)", "Status": ""})
        rows.append({"Parameter": "Reaction Energy (Delta E_rxn)", "Value": f"{neb.delta_e_reaction_kcal_mol:.2f} kcal/mol ({neb.delta_e_reaction_ev:.3f} eV)", "Status": ""})
        rows.append({"Parameter": "NEB Image Discretization", "Value": f"{neb.n_images} images (Path length = {neb.total_path_length_ang:.2f} A)", "Status": ""})

    if report.ts_frequency:
        ts = report.ts_frequency
        freq_str = f"{ts.imaginary_frequency_cm1:.1f} cm^-1" if ts.imaginary_frequency_cm1 else "None"
        rows.append({"Parameter": "Transition State Imaginary Frequency", "Value": f"{freq_str} ({ts.n_imaginary_frequencies} imag mode)", "Status": ts.status})
        if ts.irc_confirmed is not None:
            rows.append({"Parameter": "Intrinsic Reaction Coordinate (IRC)", "Value": f"{'Reported by the user as connecting the intended minima (not checked)' if ts.irc_confirmed else 'Reported as not connecting the intended minima'}", "Status": "NOT_CHECKED" if ts.irc_confirmed else "WARNING"})

    if report.tst_kinetics:
        tst = report.tst_kinetics
        rows.append({"Parameter": "TST Rate Constant (298.15 K)", "Value": f"k = {tst.k_298_s_minus_1:.2e} s^-1 (t_1/2 = {tst.half_life_298_s:.2e} s)", "Status": tst.status})
        rows.append({"Parameter": "Arrhenius Parameters (A / E_a)", "Value": f"A = {tst.arrhenius_pre_exponential_a_s_minus_1:.2e} s^-1, E_a = {tst.arrhenius_e_activation_kcal_mol:.2f} kcal/mol", "Status": ""})

    if report.tunneling:
        tun = report.tunneling
        rows.append({"Parameter": "Quantum Tunneling (298.15 K)", "Value": f"kappa_Eckart = {tun.kappa_eckart_298:.2f} (kappa_Wigner = {tun.kappa_wigner_298:.2f})", "Status": tun.status})

    df_summary = pd.DataFrame(rows)

    # CSV
    csv_path = os.path.join(output_dir, "nebcert_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path

    # LaTeX
    tex_path = os.path.join(output_dir, "nebcert_summary_table.tex")
    tex_content = df_summary.to_latex(index=False, escape=False)
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("% NEBCert summary table\n")
        f.write(tex_content)
    generated["summary_tex"] = tex_path

    # 2. Methods Text: states what was checked and what was found, including warnings and failures.
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    rxn_str = meta.get("reaction", "the reaction")
    soft_str = meta.get("software", "the electronic-structure program")
    func_str = meta.get("functional", "an unstated level of theory")

    parts = [f"The reaction path of {rxn_str} was computed with {soft_str} at {func_str} "
             f"and checked with NEBCert v{__version__}."]
    if report.neb_profile:
        neb = report.neb_profile
        src = {"optimised TS": "the optimised transition state",
               "climbing image": "the converged climbing image"}.get(
            neb.barrier_source, "the maximum of a cubic spline through the band")
        conv = {True: "converged", False: "did not converge", None: "has no reported convergence status"}[neb.band_converged]
        parts.append(f"The nudged elastic band ({neb.n_images} images) {conv}; the forward barrier, taken from {src}, "
                     f"is {neb.e_forward_barrier_kcal_mol:.2f} kcal/mol and the reaction energy "
                     f"{neb.delta_e_reaction_kcal_mol:.2f} kcal/mol (electronic energies). Band check: {neb.status}.")
    if report.ts_frequency:
        ts = report.ts_frequency
        irc = {True: " An IRC connecting the intended minima was reported (not checked by NEBCert).",
               False: " The IRC did not connect the intended minima.",
               None: " No IRC was performed or reported."}[ts.irc_confirmed]
        parts.append(f"The frequency calculation gave {ts.n_imaginary_frequencies} imaginary mode(s)"
                     + (f", the largest {ts.imaginary_frequency_cm1:.1f} cm^-1" if ts.imaginary_frequency_cm1 else "")
                     + f" (check: {ts.status}).{irc}")
    if report.tst_kinetics:
        tst = report.tst_kinetics
        kind = "an electronic barrier, so it is not a transition-state-theory rate constant" if tst.status != "PASS" \
            else "the Gibbs energy of activation"
        parts.append(f"The Eyring expression gives k(298.15 K) = {tst.k_298_s_minus_1:.2e} s^-1 from {kind}.")
    if report.tunneling:
        tun = report.tunneling
        parts.append(f"One-dimensional tunnelling factors at 298.15 K are kappa_Wigner = {tun.kappa_wigner_298:.2f} and "
                     f"kappa_Eckart = {tun.kappa_eckart_298:.2f} (check: {tun.status}).")
    parts.append(f"Overall status: {report.overall_status} ({report.validation_score.lower()}).")
    if report.recommendations:
        parts.append("Issues: " + " ".join(report.recommendations))
    full_methods = " ".join(parts)

    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(full_methods + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX
    bib_path = os.path.join(output_dir, "citation.bib")
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(BIBTEX)
    generated["citation_bib"] = bib_path

    return generated

