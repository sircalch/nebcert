"""
Evaluates one run with an old NEBCert checkout (run with PYTHONPATH pointing at it), following the
1.0.0 command-line path: band from a table of image energies (1.0.0 cannot read ORCA outputs), frequencies
given on the command line, Eyring rate and tunnelling from the spline barriers. Called by benchmark_neb.py.

    python legacy_eval.py '<json: energies_ev, s, forces, frequencies>'

Prints a JSON object with the overall status, the score string and the tunnelling factors at 298.15 K.
"""
import json
import sys

from nebcert.core.neb_profile import calculate_neb_profile_analysis
from nebcert.core.ts_frequency import verify_ts_frequency_and_irc
from nebcert.core.tst_kinetics import calculate_eyring_tst_rates
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections
from nebcert.core.scoring import assess_reaction_pathway_quality


def main():
    d = json.loads(open(sys.argv[1]).read())
    neb = calculate_neb_profile_analysis(d["energies_ev"], d.get("s"), d.get("forces"))
    ts = verify_ts_frequency_and_irc(d["frequencies"]) if d.get("frequencies") else None
    tst = calculate_eyring_tst_rates(neb.e_forward_barrier_ev)
    tun = None
    if ts is not None and ts.imaginary_frequency_cm1:
        tun = calculate_quantum_tunneling_corrections(abs(ts.imaginary_frequency_cm1), neb.e_forward_barrier_ev,
                                                      neb.e_reverse_barrier_ev)
    rep = assess_reaction_pathway_quality({}, neb, ts, tst, tun)
    print(json.dumps({"status": rep.overall_status, "score": rep.validation_score,
                      "kappa_eckart": tun.kappa_eckart_298 if tun else None,
                      "kappa_wigner": tun.kappa_wigner_298 if tun else None,
                      "messages": rep.recommendations}))


if __name__ == "__main__":
    main()
