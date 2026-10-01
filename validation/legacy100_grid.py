"""
Tunnelling factors of NEBCert 1.0.0 (the only release published before 1.2.0) on the grid of
compare_overreact.py. Run with PYTHONPATH pointing at a checkout of tag v1.0.0:

    PYTHONPATH=<v1.0.0 checkout> python validation/legacy100_grid.py

Output: validation/results/eckart_v100_grid.csv (same rows as overreact_grid.csv).
"""
import os

import pandas as pd

from nebcert import __version__
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL_PER_EV = 23.060547830619


def main():
    assert __version__ == "1.0.0", __version__
    g = pd.read_csv(os.path.join(HERE, "results", "overreact_grid.csv"))[["nu", "v1", "v2", "T"]]
    out = []
    for (nu, v1, v2), grp in g.groupby(["nu", "v1", "v2"], sort=False):
        r = calculate_quantum_tunneling_corrections(nu, v1 / KCAL_PER_EV, v2 / KCAL_PER_EV,
                                                    temperatures_k=list(grp["T"]))
        for p in r.temperature_points:
            out.append({"nu": nu, "v1": v1, "v2": v2, "T": p.temperature_k,
                        "kappa_eckart_v100": p.kappa_eckart, "kappa_wigner_v100": p.kappa_wigner})
    pd.DataFrame(out).to_csv(os.path.join(HERE, "results", "eckart_v100_grid.csv"), index=False)
    print(len(out), "points")


if __name__ == "__main__":
    main()
