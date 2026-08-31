"""
Parsers for Gaussian 16 IRC calculations.
"""

from typing import Dict, Any, List, Optional
import os
import re
import numpy as np

HARTREE_TO_EV = 27.211386245988


def parse_gaussian_irc_output(filepath: str) -> Dict[str, Any]:
    """
    Parses Gaussian 16 IRC calculation log.
    Extracts reaction coordinates, energies along the path, and verifies convergence.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    data : dict
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    points = []
    is_converged = False

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # Summary of reaction path:
            # Point Number:   1  Path Number:   1  Net Reaction Coordinate: -1.23456
            # Energy = -123.456789
            if "Net Reaction Coordinate =" in line:
                m_rc = re.search(r"Coordinate\s*=\s*([-\d\.]+)", line)
                if m_rc:
                    rc = float(m_rc.group(1))
                    points.append({"rc": rc, "energy_eh": None})
            elif "Energy =" in line and points and points[-1]["energy_eh"] is None:
                m_e = re.search(r"Energy\s*=\s*([-\d\.]+)", line)
                if m_e:
                    points[-1]["energy_eh"] = float(m_e.group(1))
            elif "Reaction path following is complete" in line or "Optimization completed" in line:
                is_converged = True

    valid_pts = [p for p in points if p["energy_eh"] is not None]
    valid_pts = sorted(valid_pts, key=lambda x: x["rc"])

    coordinates = [p["rc"] for p in valid_pts]
    energies_ev = [float(p["energy_eh"] * HARTREE_TO_EV) for p in valid_pts]

    return {
        "coordinates": coordinates,
        "energies_ev": energies_ev,
        "is_converged": is_converged,
        "n_points": len(valid_pts)
    }
