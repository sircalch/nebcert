"""
Quickstart API tutorial for NEBCert.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nebcert import (
    calculate_neb_profile_analysis,
    verify_ts_frequency_and_irc,
    calculate_eyring_tst_rates,
    calculate_quantum_tunneling_corrections,
    assess_reaction_pathway_quality
)
from nebcert.parsers import parse_mep_csv
from nebcert.reporters import (
    generate_nebcert_figures,
    generate_nebcert_manuscript_assets,
    generate_nebcert_html_report
)
from generate_sample_neb_data import generate_sample_neb_data


def main():
    print("Running NEBCert Python API quickstart tutorial...")
    raw_dir = "sample_neb_dataset"
    generate_sample_neb_data(raw_dir)

    out_dir = "quickstart_nebcert_output"
    os.makedirs(out_dir, exist_ok=True)

    # 1. Parse MEP CSV
    m_data = parse_mep_csv(os.path.join(raw_dir, "h_abstraction_neb.csv"))
    neb_res = calculate_neb_profile_analysis(
        energies_ev=m_data["energies_ev"],
        coordinates_s_ang=m_data["coordinates_s_ang"],
        tangent_forces_ev_ang=m_data["tangent_forces"]
    )

    # 2. Verify TS frequencies
    ts_res = verify_ts_frequency_and_irc(
        frequencies_cm1=[-1250.0, 120.0, 300.0, 1500.0],
        irc_confirmed=True
    )

    # 3. Calculate kinetics & tunneling
    tst_res = calculate_eyring_tst_rates(e_activation_ev=neb_res.e_forward_barrier_ev)
    tun_res = calculate_quantum_tunneling_corrections(
        imaginary_freq_cm1=1250.0,
        e_forward_barrier_ev=neb_res.e_forward_barrier_ev,
        e_reverse_barrier_ev=neb_res.e_reverse_barrier_ev
    )

    # 4. Consolidate report
    report = assess_reaction_pathway_quality(
        metadata={"reaction": "CH4 + OH -> CH3 + H2O", "functional": "wB97X-D3", "software": "ORCA"},
        neb_res=neb_res,
        ts_freq_res=ts_res,
        tst_res=tst_res,
        tunneling_res=tun_res
    )

    print(f"\nOverall Pathway Certification: {report.overall_status}")
    print(f"Forward Barrier: {report.neb_profile.e_forward_barrier_kcal_mol:.2f} kcal/mol ({report.neb_profile.e_forward_barrier_ev:.3f} eV)")
    print(f"TST Rate (298 K): {report.tst_kinetics.k_298_s_minus_1:.2e} s^-1")
    print(f"Eckart Tunneling Factor: {report.tunneling.kappa_eckart_298:.2f}")

    # 5. Export deliverables
    generate_nebcert_figures(report, out_dir)
    assets = generate_nebcert_manuscript_assets(report, out_dir)
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(out_dir, "report.html")
    generate_nebcert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print(f"\nCompleted! HTML report available at: {os.path.abspath(html_p)}")


if __name__ == "__main__":
    main()
