"""
High-precision reference for the Eckart tunnelling factor, used to settle differences between NEBCert and
overreact (compare_overreact.py). The transmission probability of Johnston & Heicklen (1962) is coded
independently of NEBCert and integrated with mpmath at 40 significant digits:

    kappa(T) = exp(V1/kT)/kT * integral_{max(0, V1-V2)}^{inf} P(E) exp(-E/kT) dE

    pip install mpmath
    python validation/eckart_highprecision.py

Output: validation/results/eckart_highprecision.csv (all 720 grid points).
"""
import os

import mpmath as mp
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 40
H, C, KB, EV, KCAL_PER_EV = 6.62607015e-34, 2.99792458e10, 1.380649e-23, 1.602176634e-19, 23.060547830619


def transmission(e, v1, v2, nu):
    a_ = v1 - v2
    b_ = (mp.sqrt(v1) + mp.sqrt(v2)) ** 2
    hnu = mp.mpf(H * C * nu / EV)
    c_ = hnu ** 2 * b_ / (16 * v1 * v2)
    if e <= 0 or e <= a_:
        return mp.mpf(0)
    a = mp.pi * mp.sqrt(e / c_)
    b = mp.pi * mp.sqrt((e - a_) / c_)
    d = (b_ - c_) / c_
    big = mp.cosh(mp.pi * mp.sqrt(d)) if d >= 0 else mp.cos(mp.pi * mp.sqrt(-d))
    return (mp.cosh(a + b) - mp.cosh(a - b)) / (mp.cosh(a + b) + big)


def kappa(nu, v1_kcal, v2_kcal, t):
    v1, v2 = mp.mpf(v1_kcal) / KCAL_PER_EV, mp.mpf(v2_kcal) / KCAL_PER_EV
    kt = mp.mpf(KB * t / EV)
    lo = max(mp.mpf(0), v1 - v2)
    f = lambda e: transmission(e, v1, v2, nu) * mp.exp(-(e - v1) / kt) / kt  # noqa: E731
    pts = sorted({lo, (lo + v1) / 2, v1, v1 + kt, v1 + 5 * kt, v1 + 80 * kt})
    return mp.quad(f, pts)


def main():
    d = pd.read_csv(os.path.join(HERE, "results", "overreact_grid.csv"))
    d["kappa_eckart_reference"] = [float(kappa(r.nu, r.v1, r.v2, r.T)) for r in d.itertuples()]
    d["rel_nebcert_ref"] = d.kappa_eckart_nebcert / d.kappa_eckart_reference - 1
    d["rel_overreact_ref"] = d.kappa_eckart_overreact / d.kappa_eckart_reference - 1
    d.to_csv(os.path.join(HERE, "results", "eckart_highprecision.csv"), index=False)
    print(f"{len(d)} points; max |rel| vs reference: NEBCert {d.rel_nebcert_ref.abs().max():.2e}, "
          f"overreact {d.rel_overreact_ref.abs().max():.2e}; overreact > 1e-3 at {(d.rel_overreact_ref.abs() > 1e-3).sum()} points")


if __name__ == "__main__":
    main()
