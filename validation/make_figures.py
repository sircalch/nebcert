"""
Figures and tables for the NEBCert v2 manuscript, built only from validation/results/ and the ORCA outputs.

    python validation/make_figures.py
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

RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
MM = 1 / 25.4
DOUBLE = 174 * MM
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
NEW, OLD, REF, AUX = "#2a78d6", "#e87ba4", "#1baf7a", "#eda100"

LABEL = {"h_hcl": "H + HCl", "oh_h2": "OH + H$_2$", "ch3_h2": "CH$_3$ + H$_2$",
         "sn2_f_ch3cl": "F$^-$ + CH$_3$Cl", "sn2_cl_ch3cl": "Cl$^-$ + CH$_3$Cl", "hcoh_to_h2co": "HCOH $\\to$ H$_2$CO",
         "h2co_dissociation": "H$_2$CO $\\to$ H$_2$ + CO", "hcn_neb": "HCN $\\to$ HNC",
         "nh3_inversion_neb": "NH$_3$ inversion"}
TEXLABEL = {k: v.replace("$\\to$", "$\\rightarrow$") for k, v in LABEL.items()}


def setup():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2, "lines.markersize": 4,
        "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42})


def panel(ax, letter, x=-0.17):
    ax.text(x, 1.03, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def fig_eckart():
    a = pd.read_csv(os.path.join(RES, "eckart_highprecision.csv"))
    b = pd.read_csv(os.path.join(RES, "eckart_v100_grid.csv"))
    m = a.merge(b, on=["nu", "v1", "v2", "T"])
    fig, (p, q) = plt.subplots(1, 2, figsize=(DOUBLE, 68 * MM))
    floor = 1e-16
    x = m.kappa_eckart_reference
    p.scatter(x, np.maximum(m.rel_overreact_ref.abs(), floor), s=9, marker="s", facecolor="none",
              edgecolor=AUX, linewidth=0.6, label="overreact")
    p.scatter(x, np.maximum(m.rel_nebcert_ref.abs(), floor), s=7, color=NEW, linewidth=0, label="NEBCert 1.2.0")
    p.set_xscale("log")
    p.set_yscale("log")
    p.set_ylim(floor / 3, 1)
    p.set_xlabel(r"$\kappa_\mathrm{Eckart}$ (reference)")
    p.set_ylabel(r"$|\kappa/\kappa_\mathrm{ref} - 1|$")
    p.axhline(1e-3, color=INK2, lw=0.6, ls=":")
    p.legend(loc="upper right", handletextpad=0.3)
    panel(p, "a")
    r = m.kappa_eckart_v100 / m.kappa_eckart_reference
    q.scatter(x, r, s=7, color=OLD, linewidth=0, label="NEBCert 1.0.0")
    q.axhline(1, color=INK2, lw=0.8)
    q.set_xscale("log")
    q.set_yscale("log")
    q.set_xlabel(r"$\kappa_\mathrm{Eckart}$ (reference)")
    q.set_ylabel(r"$\kappa_\mathrm{1.0.0}/\kappa_\mathrm{ref}$")
    q.legend(loc="upper left", handletextpad=0.3)
    panel(q, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig2_eckart")
    return m


RUNS = ["h_hcl", "oh_h2", "ch3_h2", "sn2_cl_ch3cl", "hcoh_to_h2co", "h2co_dissociation", "hcn_neb",
        "nh3_inversion_neb"]
DEGRADED = {"h2co_dissociation_3img": "H$_2$CO $\\to$ H$_2$ + CO,\nNEB-CI with 3 images",
            "oh_h2_maxiter3": "OH + H$_2$,\nstopped after 3 iterations",
            "sn2_f_ch3cl_maxiter3": "F$^-$ + CH$_3$Cl,\nend points collapsed",
            "hcoh_to_h2co_failed": "HCOH $\\to$ H$_2$CO,\nfailed job"}


def out_path(run):
    for sub in ("neb_runs/" + run, "orca_runs", "failed_runs"):
        p = os.path.join(HERE, sub, run + ".out")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(run)


def fig_profiles(b):
    names = RUNS + list(DEGRADED)
    ncol = 4
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(DOUBLE, 44 * MM * nrow))
    for ax, run in zip(axes.flat, names):
        p = parse_orca_neb_output(out_path(run))
        row = b.loc[run]
        title = LABEL.get(run, DEGRADED.get(run, run))
        e = np.array(p["energies_eh"])
        if p["n_images_without_energy"]:
            ax.text(0.5, 0.5, f"{p['n_images_without_energy']} of {p['n_images']} images\nwithout energy",
                    transform=ax.transAxes, ha="center", va="center", fontsize=7, color=INK2)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
        else:
            s = np.array(p["coordinates_s_ang"])
            rel = (e - e[0]) * HARTREE_TO_KCAL
            prof = calculate_neb_profile_analysis(p["energies_ev"], list(s), p.get("tangent_forces"))
            if np.all(np.diff(s) > 0):
                ax.plot(prof.interpolated_s, np.array(prof.interpolated_e_rel_ev) * 23.060547830619, color=NEW, lw=1)
            ax.plot(s, rel, "o", ms=3, color=NEW, mec="white", mew=0.4)
            if p.get("ci_index") is not None and not p.get("ts") and p["neb_converged"]:
                ci = p["ci_index"]
                ax.plot(s[ci], rel[ci], "D", ms=4.5, color=AUX, mec=INK, mew=0.3)
            if p.get("ts"):
                ax.plot(s[int(np.argmax(rel))], p["ts"]["barrier_fwd_kcal"], "*", ms=7, color=AUX, mec=INK, mew=0.3)
            ax.set_xlabel("path (Å)", fontsize=6.5, labelpad=1)
            if s[-1] - s[0] < 0.05:
                ax.set_xticks([s[0], s[-1]])
                ax.set_xticklabels([f"{s[0]:.3f}", f"{s[-1]:.3f}"])
        ax.set_title(title, fontsize=7, pad=3, linespacing=1.1)
        v = {"1.2.0": row.get("verdict"), "1.1.0": row.get("verdict_v110"), "1.0.0": row.get("verdict_v100")}
        txt = "   ".join(f"{k} {short(x)}" for k, x in v.items())
        ax.text(0.0, -0.42 if not p["n_images_without_energy"] else -0.25, txt, transform=ax.transAxes,
                fontsize=5.6, color=INK2)
        ax.tick_params(labelsize=6, pad=1)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("E (kcal mol$^{-1}$)", fontsize=6.5)
    fig.tight_layout(h_pad=2.2, w_pad=1.0)
    save(fig, "fig3_profiles")


def short(v):
    if not isinstance(v, str):
        return "n/a"
    if v.startswith("ERROR"):
        return "rejected" if "no energy" in v else "crash"
    return {"WARNING": "warn", "PASS": "pass", "FAIL": "FAIL"}.get(v, v)


SRC = {"optimised TS": "", "climbing image": "$^a$", "spline": "$^b$"}


def tex_num(x):
    if x is None or pd.isna(x):
        return "--"
    if abs(x) < 1e3:
        return f"{x:.2f}"
    m, e = f"{x:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def table_runs(b):
    rows = []
    for run in RUNS + list(DEGRADED):
        r = b.loc[run]
        f = lambda k, fmt: (fmt.format(r[k]) if k in r and pd.notna(r[k]) else "--")  # noqa: E731
        rows.append(" & ".join([
            TEXLABEL.get(run, DEGRADED.get(run, run).replace("\n", " ").replace("$\\to$", "$\\rightarrow$")),
            f("barrier_spline_kcal", "{:.2f}"),
            (f("barrier_checked_kcal", "{:.2f}") + SRC.get(r.get("barrier_source"), "")) if r.get("band_status") != "FAIL" else "--",
            f("nu_imag", "{:.0f}"), tex_num(r.get("kappa_eckart_reference")), tex_num(r.get("kappa_eckart_v100")),
            short(r.get("verdict")), short(r.get("verdict_v110")), short(r.get("verdict_v100"))]) + r" \\")
    with open(os.path.join(TAB, "runs.tex"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")


def main():
    setup()
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    b = pd.read_csv(os.path.join(RES, "neb_benchmark.csv")).set_index("run")
    fig_profiles(b)
    table_runs(b)
    m = fig_eckart()
    r = m.kappa_eckart_v100 / m.kappa_eckart_reference
    print(f"Eckart grid: {len(m)} points; NEBCert max |rel| {m.rel_nebcert_ref.abs().max():.2e}; "
          f"overreact max |rel| {m.rel_overreact_ref.abs().max():.2e}; 1.0.0 ratio median {r.median():.2f}, "
          f"range {r.min():.1e}-{r.max():.1e}, >2x off at {((r > 2) | (r < 0.5)).sum()} points")


if __name__ == "__main__":
    main()
