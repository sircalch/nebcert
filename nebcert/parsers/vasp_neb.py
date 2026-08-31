"""
Parsers for VASP VTST tools: neb.dat / feppaths.
"""

from typing import Dict, Any, List, Optional
import os
import numpy as np


def parse_vasp_neb_dat(filepath: str) -> Dict[str, Any]:
    """
    Parses VASP VTST 'neb.dat' file.
    Standard columns: [image_index, distance_s_ang, relative_energy_ev, tangent_force_ev_ang, curvature].

    Parameters
    ----------
    filepath : str

    Returns
    -------
    data : dict
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    indices = []
    distances = []
    rel_energies = []
    forces = []

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_s = line.strip()
            if not line_s or line_s.startswith("#"):
                continue
            parts = line_s.split()
            if len(parts) >= 3:
                try:
                    idx = int(parts[0])
                    dist = float(parts[1])
                    e_rel = float(parts[2])
                    f_val = float(parts[3]) if len(parts) >= 4 else 0.0
                    indices.append(idx)
                    distances.append(dist)
                    rel_energies.append(e_rel)
                    forces.append(f_val)
                except ValueError:
                    continue

    if not indices:
        raise ValueError(f"No numeric NEB image data parsed from {filepath}")

    return {
        "image_indices": indices,
        "distances_ang": distances,
        "energies_ev": rel_energies,
        "tangent_forces": forces
    }
