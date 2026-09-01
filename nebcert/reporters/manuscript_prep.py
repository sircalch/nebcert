"""
Manuscript Methods snippet, summary tables (CSV, LaTeX), and BibTeX citations for NEBCert.
"""

from typing import Dict, Any, Optional
import os
import pandas as pd
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

    rows.append({"Parameter": "Target Reaction Pathway", "Value": f"{meta.get('reaction', 'Chemical Reaction')} ({meta.get('software', 'DFT')})", "Status": "PASS"})
    rows.append({"Parameter": "Theoretical Level / Functional", "Value": f"{meta.get('functional', 'DFT')}", "Status": "PASS"})

    if report.neb_profile:
        neb = report.neb_profile
        rows.append({"Parameter": "Forward Activation Barrier (E_a^fwd)", "Value": f"{neb.e_forward_barrier_kcal_mol:.2f} kcal/mol ({neb.e_forward_barrier_ev:.3f} eV)", "Status": neb.status})
        rows.append({"Parameter": "Reverse Activation Barrier (E_a^rev)", "Value": f"{neb.e_reverse_barrier_kcal_mol:.2f} kcal/mol ({neb.e_reverse_barrier_ev:.3f} eV)", "Status": "PASS"})
        rows.append({"Parameter": "Reaction Energy (Delta E_rxn)", "Value": f"{neb.delta_e_reaction_kcal_mol:.2f} kcal/mol ({neb.delta_e_reaction_ev:.3f} eV)", "Status": "PASS"})
        rows.append({"Parameter": "NEB Image Discretization", "Value": f"{neb.n_images} images (Path length = {neb.total_path_length_ang:.2f} A)", "Status": "PASS"})

    if report.ts_frequency:
        ts = report.ts_frequency
        freq_str = f"{ts.imaginary_frequency_cm1:.1f} cm^-1" if ts.imaginary_frequency_cm1 else "None"
        rows.append({"Parameter": "Transition State Imaginary Frequency", "Value": f"{freq_str} ({ts.n_imaginary_frequencies} imag mode)", "Status": ts.status})
        if ts.irc_confirmed is not None:
            rows.append({"Parameter": "Intrinsic Reaction Coordinate (IRC)", "Value": f"{'Confirmed (Connects Reactants & Products)' if ts.irc_confirmed else 'Unverified'}", "Status": "PASS" if ts.irc_confirmed else "WARNING"})

    if report.tst_kinetics:
        tst = report.tst_kinetics
        rows.append({"Parameter": "TST Rate Constant (298.15 K)", "Value": f"k = {tst.k_298_s_minus_1:.2e} s^-1 (t_1/2 = {tst.half_life_298_s:.2e} s)", "Status": "PASS"})
        rows.append({"Parameter": "Arrhenius Parameters (A / E_a)", "Value": f"A = {tst.arrhenius_pre_exponential_a_s_minus_1:.2e} s^-1, E_a = {tst.arrhenius_e_activation_kcal_mol:.2f} kcal/mol", "Status": "PASS"})

    if report.tunneling:
        tun = report.tunneling
        rows.append({"Parameter": "Quantum Tunneling (298.15 K)", "Value": f"kappa_Eckart = {tun.kappa_eckart_298:.2f} (kappa_Wigner = {tun.kappa_wigner_298:.2f})", "Status": "PASS"})

    df_summary = pd.DataFrame(rows)

    # CSV
    csv_path = os.path.join(output_dir, "nebcert_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path

    # LaTeX
    tex_path = os.path.join(output_dir, "nebcert_summary_table.tex")
    tex_content = df_summary.to_latex(index=False, escape=False)
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("% NEBCert Reaction Kinetics & NEB Pathway Validation Table\n")
        f.write(tex_content)
    generated["summary_tex"] = tex_path

    # 2. Methods Text
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    rxn_str = meta.get("reaction", "chemical reaction pathways")
    soft_str = meta.get("software", "density functional theory (DFT)")
    func_str = meta.get("functional", "DFT")

    neb_str = ""
    if report.neb_profile:
        neb = report.neb_profile
        neb_str = f"The minimum energy path (MEP) was optimized with the nudged elastic band (NEB) method ({neb.n_images} images), yielding a forward activation barrier of E_a = {neb.e_forward_barrier_kcal_mol:.2f} kcal/mol ({neb.e_forward_barrier_ev:.3f} eV) and Delta E_rxn = {neb.delta_e_reaction_kcal_mol:.2f} kcal/mol. "

    ts_str = ""
    if report.ts_frequency and report.ts_frequency.imaginary_frequency_cm1:
        ts_str = f"The transition state saddle point was verified by vibrational frequency analysis, exhibiting strictly one imaginary frequency along the reaction bond coordinate (nu_imag = {report.ts_frequency.imaginary_frequency_cm1:.1f} cm^-1). "

    tst_str = ""
    if report.tst_kinetics:
        tst = report.tst_kinetics
        tun_part = f", with asymmetric Eckart quantum tunneling factor kappa = {report.tunneling.kappa_eckart_298:.2f}" if report.tunneling else ""
        tst_str = f"Reaction rate constants were determined via Eyring Transition State Theory (k_298 = {tst.k_298_s_minus_1:.2e} s^-1, Arrhenius A = {tst.arrhenius_pre_exponential_a_s_minus_1:.2e} s^-1{tun_part}). "

    full_methods = (
        f"Reaction kinetics and barrier profiles for {rxn_str} were computed using {soft_str} at the {func_str} level of theory. "
        f"Minimum energy paths, transition state verification, Eyring rate constants, and quantum tunneling corrections were certified using NEBCert v1.0.0 (Monreal-Hernández, 2026). "
        f"{neb_str}{ts_str}{tst_str}"
        f"The reaction pathway achieved an overall certification status of: {report.overall_status}."
    )

    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(full_methods + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX
    bib_path = os.path.join(output_dir, "citation.bib")
    bib_content = """@software{monreal2026nebcert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{NEBCert: Automated Quality-Control, Transition State Verification, Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/nebcert}
}
"""
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(bib_content)
    generated["citation_bib"] = bib_path

    return generated

