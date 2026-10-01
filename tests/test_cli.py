"""
End-to-end command-line test on a real ORCA 6.1.1 NEB-TS output.
"""
import os
import subprocess
import sys

DATA = os.path.join(os.path.dirname(__file__), "data")


def run(*args, tmp):
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    return subprocess.run([sys.executable, "-m", "nebcert.cli", "assess", *args, "-o", str(tmp)],
                          capture_output=True, text=True, encoding="utf-8", env=env)


def test_cli_without_irc_flag_does_not_claim_a_failed_irc(tmp_path):
    r = run("-i", os.path.join(DATA, "hcn_neb.out"), tmp=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "Overall status: PASS" in r.stdout
    methods = open(tmp_path / "methods_snippet.txt", encoding="utf-8").read()
    assert "No IRC was performed or reported" in methods
    assert "converged climbing image" not in methods and "optimised transition state" in methods


def test_cli_rejects_table_without_energies(tmp_path):
    r = run("-i", os.path.join(DATA, "orca_neb_truncated_excerpt.out"), tmp=tmp_path)
    assert r.returncode == 1
    assert "have no energy" in r.stderr
