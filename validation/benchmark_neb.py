"""
NEBCert on real ORCA 6.1.1 NEB / NEB-CI / NEB-TS calculations, against independent references.

    python validation/benchmark_neb.py [--legacy PATH_TO_1.1.0_CHECKOUT] [--legacy100 PATH_TO_1.0.0_CHECKOUT]

Runs: validation/neb_runs/*/ (make_neb_inputs.py) and validation/orca_runs/ (HCN -> HNC, NH3 inversion).
For every run the values NEBCert reads from the .out file are compared with
  * the band energies ORCA writes at full precision to <job>.final.interp (the .out table has 5 decimals);
  * the transition-state energy ORCA writes to <job>_NEB-TS_converged.xyz;
  * the harmonic frequencies ORCA writes to <job>.hess (cclib cannot read NEB outputs);
  * the Eckart and Wigner factors of overreact and of the mpmath reference (eckart_highprecision.py),
    for the barrier heights and imaginary frequency of the run.
The command-line verdict of NEBCert (this version and, with --legacy, the unreleased 1.1.0) is recorded for every
run, including the two deliberately degraded ones. Version 1.0.0, the only one published before this release,
cannot read ORCA outputs; with --legacy100 it is given the same band and frequencies as a table
(legacy_eval.py).

Must run in an environment with overreact and mpmath. Output: validation/results/neb_benchmark.csv
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from nebcert.parsers.orca_neb import parse_orca_neb_output, HARTREE_TO_KCAL  # noqa: E402
from nebcert.core.neb_profile import calculate_neb_profile_analysis  # noqa: E402
from nebcert.core.ts_frequency import verify_ts_frequency_and_irc  # noqa: E402
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections  # noqa: E402
from eckart_highprecision import kappa as kappa_reference  # noqa: E402

from overreact import tunnel  # noqa: E402

EV_PER_HARTREE = 27.211386245988
J_PER_KCAL = 4184.0


def runs():
    out = sorted(glob.glob(os.path.join(HERE, "neb_runs", "*", "*.out")))
    out += [os.path.join(HERE, "orca_runs", f) for f in ("hcn_neb.out", "nh3_inversion_neb.out")]
    out += sorted(glob.glob(os.path.join(HERE, "failed_runs", "*.out")))
    return [o for o in out if os.path.exists(o)]


def final_interp(out):
    """Band energies (Eh, relative to image 0) of the last iteration in <job>.final.interp."""
    path = out[:-4] + ".final.interp"
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        blocks = fh.read().split("Iteration:")
    rows = []
    for ln in blocks[-1].splitlines()[2:]:
        p = ln.split()
        if len(p) != 3:
            break
        rows.append(float(p[2]))
    return rows


def ts_xyz_energy(out):
    path = out[:-4] + "_NEB-TS_converged.xyz"
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        fh.readline()
        m = re.search(r"E (-?\d+\.\d+)", fh.readline())
    return float(m.group(1)) if m else None


def hess_frequencies(out):
    path = out[:-4] + ".hess"
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = fh.read().split("$vibrational_frequencies")[1].splitlines()
    n = int(lines[1])
    f = [float(lines[2 + i].split()[1]) for i in range(n)]
    return [x for x in f if abs(x) > 1e-6]


def cli_verdict(out, pythonpath):
    env = dict(os.environ, PYTHONPATH=pythonpath, PYTHONIOENCODING="utf-8")
    tmp = os.path.join(HERE, "results", "_cli_tmp")
    r = subprocess.run([sys.executable, "-m", "nebcert.cli", "assess", "-i", out, "-o", tmp],
                       capture_output=True, text=True, env=env, cwd=pythonpath, encoding="utf-8", errors="replace")
    m = re.search(r"\[RESULT\] Overall [^:]*:\s*(\w+)", r.stdout)
    s = re.search(r"\[SCORE\]\s*(.+)", r.stdout)
    return (m.group(1) if m else "ERROR: " + (r.stderr.strip().splitlines() or ["?"])[-1]), (s.group(1).strip() if s else "")


def legacy100(p, f, checkout):
    import json
    import tempfile
    payload = {"energies_ev": p["energies_ev"], "s": p.get("coordinates_s_ang"), "forces": p.get("tangent_forces"),
               "frequencies": f}
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(payload, fh)
    env = dict(os.environ, PYTHONPATH=checkout, PYTHONIOENCODING="utf-8")
    r = subprocess.run([sys.executable, os.path.join(HERE, "legacy_eval.py"), fh.name], capture_output=True,
                       text=True, env=env, cwd=checkout)
    os.unlink(fh.name)
    if r.returncode != 0:
        last = (r.stderr.strip().splitlines() or ["?"])[-1]
        return {"status": "ERROR: " + last, "score": "", "kappa_eckart": None}
    return json.loads(r.stdout.strip().splitlines()[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy", default=None, help="checkout of NEBCert 1.1.0 (git worktree of 9903f1c)")
    ap.add_argument("--legacy100", default=None, help="checkout of NEBCert 1.0.0 (git worktree of tag v1.0.0)")
    args = ap.parse_args()
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)

    rows = []
    for out in runs():
        name = os.path.basename(out)[:-4]
        p = parse_orca_neb_output(out)
        r = {"run": name, "n_images": p["n_images"], "neb_converged": p["neb_converged"],
             "terminated_normally": "ORCA TERMINATED NORMALLY" in open(out, errors="ignore").read()[-3000:]}
        if p["energies_eh"]:
            e = np.array(p["energies_eh"])
            ref = final_interp(out)
            if ref is not None and len(ref) == len(e):
                r["band_max_abs_dev_kcal"] = float(np.max(np.abs((e - e[0]) - np.array(ref))) * HARTREE_TO_KCAL)
            prof = calculate_neb_profile_analysis(p["energies_ev"], p.get("coordinates_s_ang"), p.get("tangent_forces"))
            r["barrier_spline_kcal"] = prof.e_forward_barrier_kcal_mol
            r["barrier_hei_kcal"] = float((e.max() - e[0]) * HARTREE_TO_KCAL)
            r["delta_e_kcal"] = prof.delta_e_reaction_kcal_mol
            r["well_kcal"] = prof.deepest_intermediate_well_kcal_mol
            ts0 = p.get("ts")
            chk = calculate_neb_profile_analysis(p["energies_ev"], p.get("coordinates_s_ang"), p.get("tangent_forces"),
                                                 band_converged=p["neb_converged"],
                                                 ts_energy_rel_ev=ts0["barrier_fwd_ev"] if ts0 else None,
                                                 ts_optimisation_converged=ts0["converged"] if ts0 else None,
                                                 climbing_image_index=p["ci_index"])
            r["barrier_checked_kcal"] = chk.e_forward_barrier_kcal_mol
            r["barrier_source"] = chk.barrier_source
            r["band_status"] = chk.status
            r["band_message"] = chk.diagnostic_message
        ts = p.get("ts")
        if ts:
            r["ts_converged"] = ts["converged"]
            r["barrier_ts_kcal"] = ts["barrier_fwd_kcal"]
            r["reverse_ts_kcal"] = ts["barrier_rev_ev"] * 23.060547830619
            e_xyz = ts_xyz_energy(out)
            if e_xyz is not None:
                r["ts_energy_dev_kcal"] = (ts["energy_eh"] - e_xyz) * HARTREE_TO_KCAL
        f = p.get("frequencies")
        if f:
            ff = verify_ts_frequency_and_irc(f)
            r["n_imag"] = ff.n_imaginary_frequencies
            r["nu_imag"] = ff.imaginary_frequency_cm1
            r["ts_freq_status"] = ff.status
            fref = hess_frequencies(out)
            if fref is not None and len(fref) == len(f):
                r["freq_max_abs_dev_hess"] = float(np.max(np.abs(np.array(f) - np.array(fref))))
        if ts and f and r.get("nu_imag") and r["barrier_ts_kcal"] > 0 and r["reverse_ts_kcal"] > 0:
            nu, v1, v2 = abs(r["nu_imag"]), r["barrier_ts_kcal"], r["reverse_ts_kcal"]
            tun = calculate_quantum_tunneling_corrections(nu, v1 / 23.060547830619, v2 / 23.060547830619,
                                                          temperatures_k=[298.15])
            r["kappa_eckart_nebcert"] = tun.kappa_eckart_298
            r["kappa_eckart_overreact"] = float(tunnel.eckart(nu, v1 * J_PER_KCAL, v2 * J_PER_KCAL, temperature=298.15))
            r["kappa_eckart_reference"] = float(kappa_reference(nu, v1, v2, 298.15))
            r["kappa_wigner_nebcert"] = tun.kappa_wigner_298
            r["kappa_wigner_overreact"] = float(tunnel.wigner(nu, temperature=298.15))
            r["tunnel_status"] = tun.status
        r["verdict"], r["score"] = cli_verdict(out, ROOT)
        if args.legacy:
            r["verdict_v110"], r["score_v110"] = cli_verdict(out, os.path.abspath(args.legacy))
        if args.legacy100 and p["energies_eh"]:
            l1 = legacy100(p, f, os.path.abspath(args.legacy100))
            r["verdict_v100"], r["score_v100"] = l1["status"], l1["score"]
            r["kappa_eckart_v100"] = l1["kappa_eckart"]
        rows.append(r)
        print(name, r.get("verdict"), r.get("verdict_v110", ""), r.get("verdict_v100", ""), flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(HERE, "results", "neb_benchmark.csv"), index=False)
    pd.set_option("display.width", 250)
    print(d.drop(columns=[c for c in d.columns if c.startswith("score")]).T.to_string())


if __name__ == "__main__":
    main()
