"""
Graphical abstract (square, 5 cm) for the NEBCert v2 manuscript: the HCN -> HNC band from ORCA with the
Eckart factor of versions 1.0.0 and 1.2.0 against the 40-digit reference.

    python validation/make_graphical_abstract.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from nebcert.parsers.orca_neb import parse_orca_neb_output, HARTREE_TO_KCAL  # noqa: E402
from nebcert.core.neb_profile import calculate_neb_profile_analysis  # noqa: E402

INK, INK2, NEW, OLD, AUX = "#0b0b0b", "#52514e", "#2a78d6", "#e87ba4", "#eda100"


def main():
    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                                "font.size": 6.5, "axes.spines.top": False, "axes.spines.right": False,
                                "axes.linewidth": 0.6, "pdf.fonttype": 42})
    b = pd.read_csv(os.path.join(HERE, "results", "neb_benchmark.csv")).set_index("run").loc["hcn_neb"]
    p = parse_orca_neb_output(os.path.join(HERE, "orca_runs", "hcn_neb.out"))
    s = np.array(p["coordinates_s_ang"])
    e = (np.array(p["energies_eh"]) - p["energies_eh"][0]) * HARTREE_TO_KCAL
    prof = calculate_neb_profile_analysis(p["energies_ev"], list(s), p["tangent_forces"])

    fig, ax = plt.subplots(figsize=(5 / 2.54 * 2, 5 / 2.54 * 2))  # drawn at 2x, i.e. 10 cm square
    ax.plot(prof.interpolated_s, np.array(prof.interpolated_e_rel_ev) * 23.060547830619, color=NEW, lw=1.6)
    ax.plot(s, e, "o", ms=5, color=NEW, mec="white", mew=0.6)
    ax.plot(s[int(np.argmax(e))], b.barrier_ts_kcal, "*", ms=14, color=AUX, mec=INK, mew=0.5)
    ax.set_xlabel("path (Å)", fontsize=9)
    ax.set_ylabel("E (kcal mol$^{-1}$)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_ylim(-3, 78)
    ax.set_title("HCN → HNC, ORCA NEB-TS", fontsize=10, color=INK, pad=6)
    ax.text(0.03, 0.97, "Eckart κ(298 K)", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
    ax.text(0.03, 0.88, f"reference   {b.kappa_eckart_reference:.2f}", transform=ax.transAxes, fontsize=9, va="top")
    ax.text(0.03, 0.80, f"NEBCert 1.2.0   {b.kappa_eckart_nebcert:.2f}", transform=ax.transAxes, fontsize=9,
            va="top", color=NEW)
    m, x = f"{b.kappa_eckart_v100:.1e}".split("e")
    ax.text(0.03, 0.72, f"NEBCert 1.0.0   {m}×10$^{{{int(x)}}}$", transform=ax.transAxes, fontsize=9, va="top",
            color=OLD)
    ax.text(0.97, 0.05, "failed job: 1.0.0 PASS\n1.2.0 rejected", transform=ax.transAxes, fontsize=8, ha="right",
            va="bottom", color=INK2)
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "graphical_abstract")
    fig.savefig(out + ".png", dpi=600)
    fig.savefig(out + ".pdf")
    print("written", out)


if __name__ == "__main__":
    main()
