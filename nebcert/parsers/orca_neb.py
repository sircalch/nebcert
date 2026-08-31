"""
Parsers for ORCA NEB calculations.
"""

from typing import Dict, Any, List, Optional
import os
import re
import numpy as np

HARTREE_TO_EV = 27.211386245988


def parse_orca_neb_output(filepath: str) -> Dict[str, Any]:
    """
    Parses ORCA NEB output log to extract image energies, MEP summary, and TS barrier.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    data : dict
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    energies_eh = []
    forces_max = None
    is_converged = False

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # Search for final energy summary table:
    # Image   Energy (Eh)   Distance (Bohr)
    m = re.findall(r"Image\s+(\d+)\s+:\s+([-\d\.]+)\s+Eh", content)
    if m:
        for idx_str, e_str in m:
            energies_eh.append(float(e_str))

    if "THE NEB OPTIMIZATION HAS CONVERGED" in content:
        is_converged = True

    energies_ev = [float(e * HARTREE_TO_EV) for e in energies_eh] if energies_eh else []

    return {
        "energies_ev": energies_ev,
        "n_images": len(energies_ev),
        "is_converged": is_converged
    }
