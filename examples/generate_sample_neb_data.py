"""
Generates sample NEB reaction pathway CSV and VASP VTST neb.dat files.
"""

import os
import pandas as pd
import numpy as np


def generate_sample_neb_data(output_dir: str = "sample_neb_dataset"):
    os.makedirs(output_dir, exist_ok=True)

    # 1. 9-image NEB reaction path CSV
    s_coords = [0.0, 0.35, 0.70, 1.05, 1.40, 1.75, 2.10, 2.45, 2.80]
    energies_rel_ev = [0.000, 0.045, 0.120, 0.195, 0.225, 0.110, -0.220, -0.510, -0.620]
    tangent_forces = [0.00, 0.02, 0.03, 0.02, 0.01, 0.03, 0.02, 0.02, 0.00]

    df = pd.DataFrame({
        "image": list(range(len(s_coords))),
        "distance": s_coords,
        "energy": energies_rel_ev,
        "force": tangent_forces
    })
    csv_p = os.path.join(output_dir, "h_abstraction_neb.csv")
    df.to_csv(csv_p, index=False)

    # 2. VASP neb.dat format
    dat_p = os.path.join(output_dir, "neb.dat")
    with open(dat_p, "w", encoding="utf-8") as f:
        f.write("# Image  Distance_A   Energy_eV    Force_eV_A   Curvature\n")
        for i, (s, e, f_val) in enumerate(zip(s_coords, energies_rel_ev, tangent_forces)):
            f.write(f"  {i:3d}    {s:10.5f}   {e:10.5f}   {f_val:10.5f}   0.00000\n")

    print(f"Generated sample NEB dataset at: {os.path.abspath(output_dir)}/")


if __name__ == "__main__":
    generate_sample_neb_data()
