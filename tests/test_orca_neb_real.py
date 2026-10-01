import os
import pytest
from nebcert.parsers.orca_neb import parse_orca_neb_output

F = os.path.join(os.path.dirname(__file__), "data", "nh3_inversion_neb.out")


def test_orca6_neb_ts_ammonia_inversion():
    d = parse_orca_neb_output(F)
    assert d["n_images"] == 9 and d["ci_index"] == 4 and d["neb_converged"]
    assert d["ts"]["barrier_fwd_kcal"] == pytest.approx(5.93, abs=0.02)   # direct TS - min: 5.93 kcal/mol
    assert d["ts"]["converged"] is True
    assert d["frequencies"][0] == pytest.approx(-830.15)
    assert sum(f < 0 for f in d["frequencies"]) == 1


def test_orca6_neb_ts_hcn_to_hnc():
    d = parse_orca_neb_output(os.path.join(os.path.dirname(__file__), "data", "hcn_neb.out"))
    assert d["n_images"] == 11 and d["ci_index"] == 4 and d["neb_converged"]
    assert d["ts"]["converged"] is True        # banner reads "THE TS OPTIMIZATION HAS CONVERGED"
    assert d["ts"]["barrier_fwd_kcal"] == pytest.approx(47.78, abs=0.02)
    assert d["energies_ev"][-1] * 23.0605 == pytest.approx(13.56, abs=0.05)   # HNC above HCN
    assert d["frequencies"] == pytest.approx([-1116.47, 2089.78, 2625.04])


def test_truncated_output_reports_images_without_energy():
    d = parse_orca_neb_output(os.path.join(os.path.dirname(F), "orca_neb_truncated_excerpt.out"))
    assert d["n_images"] == 10
    assert d["n_images_without_energy"] == 9
    assert d["neb_converged"] is False
