"""
Tests for scoring, manuscript reporting, and CLI demo execution in NEBCert.
"""

import os
import tempfile
import numpy as np
import pytest
from nebcert.core.neb_profile import calculate_neb_profile_analysis
from nebcert.core.ts_frequency import verify_ts_frequency_and_irc
from nebcert.core.tst_kinetics import calculate_eyring_tst_rates
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections
from nebcert.core.scoring import assess_reaction_pathway_quality
from nebcert.reporters.plot_generator import generate_nebcert_figures
from nebcert.reporters.manuscript_prep import generate_nebcert_manuscript_assets
from nebcert.reporters.html_report import generate_nebcert_html_report
from nebcert.cli import run_demo


def test_full_nebcert_validation_pipeline():
    meta = {
        "reaction": "A + B -> [TS] -> C",
        "functional": "B3LYP",
        "software": "Gaussian"
    }

    neb_res = calculate_neb_profile_analysis(
        energies_ev=[0.0, 0.2, 0.45, 0.1, -0.2],
        coordinates_s_ang=[0.0, 0.5, 1.0, 1.5, 2.0],
        tangent_forces_ev_ang=[0.0, 0.02, 0.01, 0.02, 0.0]
    )

    ts_res = verify_ts_frequency_and_irc(
        frequencies_cm1=[-950.0, 120.0, 300.0, 1500.0],
        irc_confirmed=True
    )

    tst_res = calculate_eyring_tst_rates(
        e_activation_ev=neb_res.e_forward_barrier_ev
    )

    tun_res = calculate_quantum_tunneling_corrections(
        imaginary_freq_cm1=950.0,
        e_forward_barrier_ev=neb_res.e_forward_barrier_ev,
        e_reverse_barrier_ev=neb_res.e_reverse_barrier_ev
    )

    report = assess_reaction_pathway_quality(
        metadata=meta,
        neb_res=neb_res,
        ts_freq_res=ts_res,
        tst_res=tst_res,
        tunneling_res=tun_res
    )

    assert report.overall_status == "PASS"

    with tempfile.TemporaryDirectory() as tmpdir:
        plots = generate_nebcert_figures(report, tmpdir, formats=["png", "svg"])
        assert len(plots) > 0
        for p in plots:
            assert os.path.exists(p)

        assets = generate_nebcert_manuscript_assets(report, tmpdir)
        assert os.path.exists(assets["summary_csv"])
        assert os.path.exists(assets["summary_tex"])
        assert os.path.exists(assets["methods_text"])
        assert os.path.exists(assets["citation_bib"])

        html_p = os.path.join(tmpdir, "report.html")
        generate_nebcert_html_report(report, html_p, methods_text="Sample methods", citation_bib="@software{}")
        assert os.path.exists(html_p)
        assert os.path.getsize(html_p) > 500


def test_cli_demo_execution():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_demo(output_dir=tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "report.html"))
        assert os.path.exists(os.path.join(tmpdir, "nebcert_summary_table.csv"))
        assert os.path.exists(os.path.join(tmpdir, "citation.bib"))
