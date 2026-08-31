"""
Tests for parsers (VASP, ORCA, Gaussian, CSV) in NEBCert.
"""

import os
import tempfile
import numpy as np
import pytest
from nebcert.parsers.vasp_neb import parse_vasp_neb_dat
from nebcert.parsers.generic_mep_csv import parse_mep_csv


def test_vasp_neb_dat_parser():
    content = """# Image  Distance    Energy       Force      Curvature
    0     0.00000     0.00000     0.00000    0.00000
    1     0.50000     0.25000     0.02000   -0.50000
    2     1.00000     0.60000     0.01000   -1.20000
    3     1.50000     0.30000     0.01500   -0.40000
    4     2.00000    -0.20000     0.00000    0.00000
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
        f.write(content)
        f_path = f.name

    try:
        data = parse_vasp_neb_dat(f_path)
        assert len(data["image_indices"]) == 5
        assert np.isclose(data["energies_ev"][2], 0.60)
        assert np.isclose(data["distances_ang"][4], 2.0)
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_generic_mep_csv_parser():
    content = """img,energy,dist,force
0,0.0,0.0,0.00
1,0.3,0.5,0.02
2,0.7,1.0,0.01
3,0.2,1.5,0.03
4,-0.1,2.0,0.00
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        f_path = f.name

    try:
        data = parse_mep_csv(f_path)
        assert data["n_images"] == 5
        assert len(data["energies_ev"]) == 5
        assert data["coordinates_s_ang"] == [0.0, 0.5, 1.0, 1.5, 2.0]
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)
