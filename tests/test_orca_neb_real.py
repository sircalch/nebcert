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
