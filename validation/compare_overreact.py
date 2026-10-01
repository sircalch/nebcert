"""
Tunnelling and rate constants of NEBCert against overreact (Schneider et al.; github.com/geem-lab/overreact,
commit dec053b), an independent implementation of the Eckart (Johnston & Heicklen 1962) and Wigner
corrections and of the Eyring equation.

    python validation/compare_overreact.py

Grid: imaginary frequency 200-3000 cm-1, forward barrier 2-40 kcal/mol, reverse barrier 2-60 kcal/mol,
temperature 200-1000 K. Output: validation/results/overreact_grid.csv
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections  # noqa: E402
from nebcert.core.tst_kinetics import calculate_eyring_tst_rates  # noqa: E402

from overreact import rates, tunnel  # noqa: E402

KCAL_PER_EV = 23.060547830619
J_PER_KCAL = 4184.0

FREQS = [200, 500, 1000, 1500, 2000, 3000]
V1 = [2, 5, 10, 20, 40]
V2 = [2, 5, 10, 20, 40, 60]
TEMPS = [200.0, 298.15, 500.0, 1000.0]


def main():
    rows = []
    for nu, v1, v2 in itertools.product(FREQS, V1, V2):
        r = calculate_quantum_tunneling_corrections(nu, v1 / KCAL_PER_EV, v2 / KCAL_PER_EV, temperatures_k=TEMPS)
        for p in r.temperature_points:
            t = p.temperature_k
            k_e = float(tunnel.eckart(nu, v1 * J_PER_KCAL, v2 * J_PER_KCAL, temperature=t))
            k_w = float(tunnel.wigner(nu, temperature=t))
            rows.append({"nu": nu, "v1": v1, "v2": v2, "T": t, "kappa_eckart_nebcert": p.kappa_eckart,
                         "kappa_eckart_overreact": k_e, "kappa_wigner_nebcert": p.kappa_wigner,
                         "kappa_wigner_overreact": k_w})
    d = pd.DataFrame(rows)
    d["rel_eckart"] = d.kappa_eckart_nebcert / d.kappa_eckart_overreact - 1
    d["rel_wigner"] = d.kappa_wigner_nebcert / d.kappa_wigner_overreact - 1

    e_rows = []
    for v in (5.0, 10.0, 20.0, 30.0):
        r = calculate_eyring_tst_rates(v / KCAL_PER_EV, temperatures_k=TEMPS)
        for p in r.rate_points if hasattr(r, "rate_points") else r.temperature_points:
            k_ref = float(np.ravel(rates.eyring(v * J_PER_KCAL, temperature=p.temperature_k))[0])
            e_rows.append({"barrier_kcal": v, "T": p.temperature_k, "k_nebcert": p.k_tst_s_minus_1, "k_overreact": k_ref})
    e = pd.DataFrame(e_rows)
    e["rel"] = e.k_nebcert / e.k_overreact - 1

    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    d.to_csv(os.path.join(HERE, "results", "overreact_grid.csv"), index=False)
    e.to_csv(os.path.join(HERE, "results", "overreact_eyring.csv"), index=False)
    print(f"{len(d)} Eckart/Wigner points; max |rel| Eckart {d.rel_eckart.abs().max():.2e}, Wigner {d.rel_wigner.abs().max():.2e}")
    print(f"{len(e)} Eyring points; max |rel| {e.rel.abs().max():.2e}")
    worst = d.loc[d.rel_eckart.abs().sort_values(ascending=False).index[:8]]
    print(worst[["nu", "v1", "v2", "T", "kappa_eckart_nebcert", "kappa_eckart_overreact", "rel_eckart"]].to_string())


if __name__ == "__main__":
    main()
