"""
Parser for generic tabular MEP / NEB CSV/TSV files.
"""

from typing import Dict, Any, List, Tuple
import os
import pandas as pd
import numpy as np


def parse_mep_csv(filepath: str) -> Dict[str, Any]:
    """
    Parses generic MEP/NEB CSV or TSV files.
    Recognized columns: 'image'/'index', 'energy'/'energy_ev', 'coordinate'/'s'/'dist', 'force'/'tangent_force'.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    data : dict
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    sep = r"\s+" if filepath.endswith(".dat") else ("," if filepath.endswith(".csv") else None)
    df = pd.read_csv(filepath, sep=sep, engine="python" if sep is None else None)

    col_map = {}
    for c in df.columns:
        c_low = str(c).lower().strip()
        if c_low in ["image", "index", "img", "id", "point"]:
            col_map[c] = "image"
        elif c_low in ["energy", "energy_ev", "e", "e_rel", "e_ev"]:
            col_map[c] = "energy"
        elif c_low in ["coordinate", "s", "dist", "distance", "rx_coord", "rc"]:
            col_map[c] = "coordinate"
        elif c_low in ["force", "f", "tangent_force", "f_tangent", "max_force"]:
            col_map[c] = "force"

    df = df.rename(columns=col_map)

    if "energy" not in df.columns:
        raise ValueError(f"MEP/NEB table must contain an 'energy' column. Found: {list(df.columns)}")

    return {
        "energies_ev": df["energy"].astype(float).tolist(),
        "coordinates_s_ang": df["coordinate"].astype(float).tolist() if "coordinate" in df.columns else None,
        "tangent_forces": df["force"].astype(float).tolist() if "force" in df.columns else None,
        "n_images": len(df)
    }
