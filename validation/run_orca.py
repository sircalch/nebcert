"""
Runs every input in validation/neb_runs/*/ with ORCA (serial jobs, several at a time).

    python validation/run_orca.py [n_parallel]

Jobs with an output that already ends in "ORCA TERMINATED NORMALLY" are skipped.
Set ORCA to the orca executable if it is not at C:\\ORCA_6.1.1\\orca.exe.
"""
import ctypes
import glob
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ORCA = os.environ.get("ORCA", r"C:\ORCA_6.1.1\orca.exe")


def done(out):
    if not os.path.exists(out):
        return False
    with open(out, encoding="utf-8", errors="ignore") as fh:
        return "ORCA TERMINATED NORMALLY" in fh.read()[-3000:]


def run(inp):
    d, name = os.path.dirname(inp), os.path.basename(inp)[:-4]
    out = os.path.join(d, name + ".out")
    if done(out):
        return name, "skipped"
    t = time.time()
    with open(out, "w") as fh:
        subprocess.run([ORCA, name + ".inp"], cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    return name, f"{'ok' if done(out) else 'ended abnormally'} in {time.time() - t:.0f} s"


def main():
    ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)  # keep the machine awake
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    inputs = sorted(glob.glob(os.path.join(HERE, "neb_runs", "*", "*.inp")))
    # small molecules first so that results arrive early
    inputs.sort(key=lambda p: sum(1 for line in open(p) if len(line.split()) == 4))
    with ThreadPoolExecutor(max_workers=n) as ex:
        for name, status in ex.map(run, inputs):
            print(name, status, flush=True)


if __name__ == "__main__":
    main()
