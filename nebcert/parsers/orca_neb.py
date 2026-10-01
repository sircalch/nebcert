"""
Parser for ORCA (5/6) NEB, NEB-CI and NEB-TS output files.

Reads the last "PATH SUMMARY" table (image, distance, energy, max(|Fp|), RMS(Fp)), the climbing
image, the "PATH SUMMARY FOR NEB-TS" transition-state row, convergence banners and the final
vibrational frequencies (NEB-TS ... Freq).
"""

from typing import Dict, Any, List, Optional
import os
import re

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.509474063
EH_PER_BOHR_TO_EV_PER_ANG = 51.42208619083232

_NUM = r"-?\d+\.\d+"
_ROW_WITH_DIST = re.compile(rf"^\s*(\d+)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})(.*)$")
_ROW_NO_DIST = re.compile(rf"^\s*(\d+|TS)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})(.*)$")


def _last_block(content: str, header_regex: str) -> Optional[str]:
    parts = re.split(header_regex, content)
    if len(parts) < 2:
        return None
    body = parts[-1]
    # table ends at the first empty line after the column header
    lines = body.splitlines()
    out, started = [], False
    for ln in lines:
        if "Image" in ln and "E(Eh)" in ln:
            started = True
            continue
        if started:
            if not ln.strip():
                break
            out.append(ln)
    return "\n".join(out)


def _final_frequencies(content: str) -> Optional[List[float]]:
    blocks = re.split(r"VIBRATIONAL FREQUENCIES\s*\n-+\s*\n", content)
    if len(blocks) < 2:
        return None
    body = re.split(r"NORMAL MODES|IR SPECTRUM", blocks[-1])[0]
    freqs = [float(v) for v in re.findall(r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1", body, re.M)]
    freqs = [f for f in freqs if abs(f) > 1e-6]
    return freqs or None


def parse_orca_neb_output(filepath: str) -> Dict[str, Any]:
    """
    Parses an ORCA NEB / NEB-CI / NEB-TS output file.

    Returns
    -------
    data : dict
        energies_eh, energies_ev (relative to image 0), coordinates_s_ang (cumulative distance),
        tangent_forces (max |F_perp| per image, eV/Angstrom), rms_forces (eV/Angstrom),
        ci_index, n_images, n_images_without_energy (images printed with E = 0, i.e. not computed), neb_converged, is_converged (alias), ts (dict or None: energy_eh,
        barrier_fwd_ev, barrier_rev_ev, barrier_fwd_kcal, converged), frequencies (final
        vibrational frequencies, translations/rotations excluded, or None).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath, "r", encoding="utf-8", errors="ignore") as fh:
        content = fh.read()

    images, dists, energies, fmax, frms = [], [], [], [], []
    ci_index = None
    block = _last_block(content, r"\n\s*PATH SUMMARY\s*\n")
    if block:
        for ln in block.splitlines():
            m = _ROW_WITH_DIST.match(ln)
            if not m:
                continue
            images.append(int(m.group(1)))
            dists.append(float(m.group(2)))
            energies.append(float(m.group(3)))
            fmax.append(float(m.group(5)) * EH_PER_BOHR_TO_EV_PER_ANG)
            frms.append(float(m.group(6)) * EH_PER_BOHR_TO_EV_PER_ANG)
            if "CI" in m.group(7):
                ci_index = len(images) - 1

    neb_converged = bool(re.search(r"THE NEB OPTIMIZATION HAS CONVERGED", content))

    ts = None
    ts_block = _last_block(content, r"\n\s*PATH SUMMARY FOR NEB-TS\s*\n")
    if ts_block:
        for ln in ts_block.splitlines():
            m = _ROW_NO_DIST.match(ln)
            if m and m.group(1) == "TS":
                e_ts = float(m.group(2))
                ts = {"energy_eh": e_ts}
                break
    if ts is not None and energies:
        ts["barrier_fwd_ev"] = (ts["energy_eh"] - energies[0]) * HARTREE_TO_EV
        ts["barrier_rev_ev"] = (ts["energy_eh"] - energies[-1]) * HARTREE_TO_EV
        ts["barrier_fwd_kcal"] = (ts["energy_eh"] - energies[0]) * HARTREE_TO_KCAL
        after_neb = content.split("PATH SUMMARY FOR NEB-TS")[0]
        ts_part = after_neb.split("THE NEB OPTIMIZATION HAS CONVERGED")[-1] if neb_converged else after_neb
        ts["converged"] = bool(re.search(r"THE (?:TS )?OPTIMIZATION HAS CONVERGED", ts_part))

    # ORCA prints the PATH SUMMARY table once before the band is optimised, with E = 0.00000 for the images
    # not yet computed; a job whose image calculations fail (or that is killed at that stage) leaves only this table.
    n_missing = sum(1 for e in energies if e == 0.0)
    energies_ev = [(e - energies[0]) * HARTREE_TO_EV for e in energies] if energies else []
    return {
        "n_images_without_energy": n_missing,
        "energies_eh": energies,
        "energies_ev": energies_ev,
        "coordinates_s_ang": dists or None,
        "tangent_forces": fmax or None,
        "rms_forces": frms or None,
        "ci_index": ci_index,
        "n_images": len(energies),
        "neb_converged": neb_converged,
        "is_converged": neb_converged,
        "ts": ts,
        "frequencies": _final_frequencies(content),
    }
