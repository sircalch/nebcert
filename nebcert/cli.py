"""
Command Line Interface (CLI) for NEBCert.
"""

import sys
import os
import argparse
import numpy as np

from nebcert import __version__
from nebcert.parsers.vasp_neb import parse_vasp_neb_dat
from nebcert.parsers.orca_neb import parse_orca_neb_output
from nebcert.parsers.gaussian_irc import parse_gaussian_irc_output
from nebcert.parsers.generic_mep_csv import parse_mep_csv

from nebcert.core.neb_profile import calculate_neb_profile_analysis
from nebcert.core.ts_frequency import verify_ts_frequency_and_irc
from nebcert.core.tst_kinetics import calculate_eyring_tst_rates
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections
from nebcert.core.scoring import assess_reaction_pathway_quality

from nebcert.reporters.plot_generator import generate_nebcert_figures
from nebcert.reporters.manuscript_prep import generate_nebcert_manuscript_assets
from nebcert.reporters.html_report import generate_nebcert_html_report


def print_banner():
    banner = rf"""
  _   _ ______ _____   _____          _   
 | \ | |  ____|  _ \ / ____|        | |  
 |  \| | |__  | |_) | |     ___ _ __| |_ 
 | . ` |  __| |  _ <| |    / _ \ '__| __|
 | |\  | |____| |_) | |___|  __/ |  | |_ 
 |_| \_|______|____/ \_____\___|_|   \__| v{__version__}

 Transition State, NEB Reaction Pathways & Quantum Tunneling Toolkit
 Monreal-Hernández et al., 2026
"""
    print(banner)


def run_demo(output_dir: str = "nebcert_demo_output"):
    """
    Executes a benchmark demonstration evaluating a hydrogen atom abstraction reaction:
    CH4 + OH* -> CH3* + H2O with 9-image NEB band, 1st-order saddle point TS (nu_imag = -1250 cm^-1),
    Eyring TST kinetics, and asymmetric Eckart quantum tunneling.
    """
    print(f"\n[NEBCert] Running demonstration benchmark on Reaction Pathway (CH4 + OH* -> CH3* + H2O)...")
    os.makedirs(output_dir, exist_ok=True)

    metadata = {
        "reaction": "CH4 + OH* -> [TS]* -> CH3* + H2O (Hydrogen Abstraction)",
        "functional": "wB97X-D3 / def2-TZVP",
        "software": "ORCA 6.0.0 (CI-NEB / NumFreq)"
    }

    # 1. 9-image NEB reaction energy profile (eV)
    # Reactant = 0.0, Peak TS (img 4) = 0.225 eV (5.19 kcal/mol), Product (img 8) = -0.620 eV (-14.30 kcal/mol)
    s_coords = [0.0, 0.35, 0.70, 1.05, 1.40, 1.75, 2.10, 2.45, 2.80]
    energies_rel_ev = [0.000, 0.045, 0.120, 0.195, 0.225, 0.110, -0.220, -0.510, -0.620]
    tangent_forces = [0.00, 0.02, 0.03, 0.02, 0.01, 0.03, 0.02, 0.02, 0.00]

    print("  -> Analyzing NEB Minimum Energy Path, fitting cubic spline, and auditing CI-NEB forces...")
    neb_res = calculate_neb_profile_analysis(
        energies_ev=energies_rel_ev,
        coordinates_s_ang=s_coords,
        tangent_forces_ev_ang=tangent_forces
    )

    # 2. Transition state imaginary frequency
    print("  -> Verifying 1st-order saddle point vibrational frequencies (nu_imag = -1250.0 cm^-1)...")
    sample_freqs = [-1250.0, 120.0, 250.0, 480.0, 850.0, 1100.0, 1450.0, 2900.0, 3100.0, 3650.0]
    ts_res = verify_ts_frequency_and_irc(
        frequencies_cm1=sample_freqs,
        min_significant_freq_cm1=50.0,
        irc_confirmed=True
    )

    # 3. Eyring TST rate constants
    print("  -> Computing Eyring-Polanyi TST rate constants k(T) and Arrhenius parameters...")
    tst_res = calculate_eyring_tst_rates(
        e_activation_ev=neb_res.e_forward_barrier_ev
    )

    # 4. Quantum tunneling corrections (Eckart & Wigner)
    print("  -> Computing asymmetric Eckart quantum tunneling transmission factor kappa(T)...")
    tun_res = calculate_quantum_tunneling_corrections(
        imaginary_freq_cm1=1250.0,
        e_forward_barrier_ev=neb_res.e_forward_barrier_ev,
        e_reverse_barrier_ev=neb_res.e_reverse_barrier_ev
    )

    report = assess_reaction_pathway_quality(
        metadata=metadata,
        neb_res=neb_res,
        ts_freq_res=ts_res,
        tst_res=tst_res,
        tunneling_res=tun_res
    )

    print("  -> Generating publication-ready vector figures (NEB MEP Profile, Arrhenius Plot, Tunneling Curve)...")
    generate_nebcert_figures(report, output_dir)

    print("  -> Drafting manuscript Methods text snippet, summary LaTeX tables, and BibTeX citations...")
    assets = generate_nebcert_manuscript_assets(report, output_dir)

    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing interactive report to {html_p}...")
    generate_nebcert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print("\n" + "="*70)
    print(f" [RESULT] Overall Reaction Pathway Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    print(f" * Reaction Target : {report.metadata['reaction']}")
    print(f" * Activation Barrier: E_a^fwd = {report.neb_profile.e_forward_barrier_kcal_mol:.2f} kcal/mol ({report.neb_profile.e_forward_barrier_ev:.3f} eV) | Delta E_rxn = {report.neb_profile.delta_e_reaction_kcal_mol:.2f} kcal/mol")
    print(f" * Transition State : {report.ts_frequency.imaginary_frequency_cm1:.1f} cm^-1 (Strict 1st-order saddle point) | Status: {report.ts_frequency.status}")
    print(f" * Eyring TST Rate  : k(298 K) = {report.tst_kinetics.k_298_s_minus_1:.2e} s^-1 (Arrhenius A = {report.tst_kinetics.arrhenius_pre_exponential_a_s_minus_1:.2e} s^-1, E_a = {report.tst_kinetics.arrhenius_e_activation_kcal_mol:.2f} kcal/mol)")
    print(f" * Quantum Tunneling: kappa_Eckart(298 K) = {report.tunneling.kappa_eckart_298:.2f} (Tunneling acceleration = x{report.tunneling.kappa_eckart_298:.1f})")
    print("="*70)
    print(f"\nAll outputs successfully saved to: {os.path.abspath(output_dir)}/")
    print(f"Open {os.path.abspath(html_p)} in your browser to inspect the full report.\n")


def run_assess(args):
    """
    Evaluates user-provided MEP/NEB files or frequencies.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    neb_res = None
    ts_res = None
    tst_res = None
    tun_res = None

    # 1. Parse NEB / MEP input
    if args.input_neb:
        print(f"\n[NEBCert] Parsing NEB data from: {args.input_neb}...")
        if args.input_neb.endswith(".dat"):
            try:
                m_data = parse_vasp_neb_dat(args.input_neb)
            except Exception:
                m_data = parse_mep_csv(args.input_neb)
        else:
            m_data = parse_mep_csv(args.input_neb)

        neb_res = calculate_neb_profile_analysis(
            energies_ev=m_data["energies_ev"],
            coordinates_s_ang=m_data.get("coordinates_s_ang"),
            tangent_forces_ev_ang=m_data.get("tangent_forces")
        )

    # 2. Parse / evaluate TS frequencies
    if args.frequencies:
        freq_list = [float(x) for x in args.frequencies.split(",")]
        ts_res = verify_ts_frequency_and_irc(
            frequencies_cm1=freq_list,
            irc_confirmed=args.irc
        )

    # 3. Calculate kinetics and tunneling
    ea_val = None
    if neb_res:
        ea_val = neb_res.e_forward_barrier_ev
    elif args.barrier_ev:
        ea_val = float(args.barrier_ev)
    elif args.barrier_kcal:
        ea_val = float(args.barrier_kcal) / 23.06054887

    if ea_val is not None:
        tst_res = calculate_eyring_tst_rates(e_activation_ev=ea_val)
        if ts_res and ts_res.imaginary_frequency_cm1:
            tun_res = calculate_quantum_tunneling_corrections(
                imaginary_freq_cm1=ts_res.imaginary_frequency_cm1,
                e_forward_barrier_ev=ea_val,
                e_reverse_barrier_ev=neb_res.e_reverse_barrier_ev if neb_res else None
            )

    meta = {
        "reaction": args.reaction or "Chemical Reaction Pathway",
        "functional": args.functional or "DFT",
        "software": args.software or "VASP / ORCA / Gaussian"
    }

    report = assess_reaction_pathway_quality(
        metadata=meta,
        neb_res=neb_res,
        ts_freq_res=ts_res,
        tst_res=tst_res,
        tunneling_res=tun_res
    )

    print("  -> Generating publication figures...")
    generate_nebcert_figures(report, output_dir)

    print("  -> Generating manuscript text, LaTeX summary table, and BibTeX citations...")
    assets = generate_nebcert_manuscript_assets(report, output_dir)

    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing HTML quality report to {html_p}...")
    generate_nebcert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print("\n" + "="*70)
    print(f" [RESULT] Overall Quality Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    if report.neb_profile:
        print(f" * Barrier (E_a^fwd): {report.neb_profile.e_forward_barrier_kcal_mol:.2f} kcal/mol ({report.neb_profile.e_forward_barrier_ev:.3f} eV)")
    if report.tst_kinetics:
        print(f" * TST Rate (298 K) : k = {report.tst_kinetics.k_298_s_minus_1:.2e} s^-1")
    print("="*70)
    print(f"\nReport ready at: {os.path.abspath(html_p)}\n")


def print_citation():
    bib = """@software{monreal2026nebcert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{NEBCert: Automated Quality-Control, Transition State Verification, Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/nebcert}
}"""
    print("\nIf you use NEBCert in your publications, please cite:\n")
    print("APA Style:")
    print("Monreal-Hernández, A. (2026). NEBCert: Automated Quality-Control, Transition State Verification, Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification (v1.0.0). Zenodo. https://github.com/sircalch/nebcert\n")
    print("BibTeX:")
    print(bib)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="nebcert",
        description="NEBCert: Transition State, NEB Reaction Pathways, Quantum Tunneling, and TST Kinetics Certification."
    )
    parser.add_argument("-v", "--version", action="version", version=f"nebcert {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Assess command
    assess_parser = subparsers.add_parser("assess", help="Assess NEB pathway, TS frequencies, and reaction kinetics")
    assess_parser.add_argument("-i", "--input-neb", default=None, help="Path to NEB file (neb.dat, or CSV/TSV table of image energies)")
    assess_parser.add_argument("--frequencies", default=None, help="Comma-separated vibrational frequencies in cm^-1 (e.g. '-1250,150,300,800')")
    assess_parser.add_argument("--barrier-ev", type=float, default=None, help="Activation energy in eV (if not using NEB profile)")
    assess_parser.add_argument("--barrier-kcal", type=float, default=None, help="Activation energy in kcal/mol")
    assess_parser.add_argument("--irc", action="store_true", help="Flag if IRC was confirmed")
    assess_parser.add_argument("-o", "--output", default="nebcert_output", help="Directory for output report (default: nebcert_output)")
    assess_parser.add_argument("--reaction", default=None, help="Reaction description (e.g. 'CH4 + OH -> CH3 + H2O')")
    assess_parser.add_argument("--functional", default=None, help="DFT functional / level of theory")
    assess_parser.add_argument("--software", default=None, help="Software code (e.g. 'VASP', 'ORCA', 'Gaussian')")

    # Demo command
    demo_parser = subparsers.add_parser("demo", help="Run benchmark demonstration (H-abstraction 9-image NEB + TS + TST + Tunneling)")
    demo_parser.add_argument("-o", "--output", default="nebcert_demo_output", help="Output directory (default: nebcert_demo_output)")

    # Cite command
    subparsers.add_parser("cite", help="Display BibTeX and APA citation details")

    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if args.command == "assess":
        print_banner()
        run_assess(args)
    elif args.command == "demo":
        print_banner()
        run_demo(args.output)
    elif args.command == "cite":
        print_banner()
        print_citation()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

