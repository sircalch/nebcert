"""
Writes the ORCA 6.1 NEB-TS inputs of the NEBCert v2 benchmark (validation/neb_runs/<name>/).

    python validation/make_neb_inputs.py

Level: B3LYP-D3(BJ)/def2-SVP, TightSCF, NEB-TS with 8 images and pre-optimised end points; the
transition state is followed by an analytic frequency calculation (Freq). Open shells are UKS.
Reactions (atom order identical in both end points):
    h2co_dissociation   H2CO -> H2 + CO
    hcoh_to_h2co        trans-HCOH -> H2CO (1,2-H shift)
    sn2_cl_ch3cl        Cl- + CH3Cl -> ClCH3 + Cl- (identity SN2, charge -1)
    h_hcl               H + HCl -> H2 + Cl (doublet)
    oh_h2               OH + H2 -> H2O + H (doublet)
    ch3_h2              CH3 + H2 -> CH4 + H (doublet)
Deliberately degraded runs:
    h2co_dissociation_3img   NEB-CI with 3 images and loose convergence (no TS optimisation)
    sn2_f_ch3cl_maxiter3     F- + CH3Cl -> FCH3 + Cl- (charge -1), NEB-TS limited to 3 NEB iterations. The
                             limit never acts: at this level the reactant end point, started 2.6 A from C,
                             slides into the product complex during pre-optimisation, both end points become
                             the same structure and ORCA reports "No barrier was found". Kept as a real
                             degenerate case. (A run of the same reaction without the limit died silently
                             during pre-optimisation and was not repeated.)
    oh_h2_maxiter3           NEB-TS stopped after 3 NEB iterations
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "neb_runs")

R = {
    "h2co_dissociation": (0, 1, [("C", 0, 0, 0), ("O", 0, 0, 1.21), ("H", 0, 0.94, -0.58), ("H", 0, -0.94, -0.58)],
                          [("C", 0, 0, 0), ("O", 0, 0, 1.13), ("H", 0, 0.37, -2.60), ("H", 0, -0.37, -2.60)]),
    "hcoh_to_h2co": (0, 1, [("C", 0, 0, 0), ("O", 1.32, 0, 0), ("H", -0.55, 0.95, 0), ("H", 1.65, -0.90, 0)],
                     [("C", 0, 0, 0), ("O", 1.21, 0, 0), ("H", -0.58, 0.94, 0), ("H", -0.58, -0.94, 0)]),
    "sn2_f_ch3cl": (-1, 1, [("F", 0, 0, -2.60), ("C", 0, 0, 0), ("Cl", 0, 0, 1.80),
                            ("H", 1.03, 0, -0.36), ("H", -0.515, 0.892, -0.36), ("H", -0.515, -0.892, -0.36)],
                    [("F", 0, 0, -1.40), ("C", 0, 0, 0), ("Cl", 0, 0, 3.10),
                     ("H", 1.03, 0, 0.36), ("H", -0.515, 0.892, 0.36), ("H", -0.515, -0.892, 0.36)]),
    "sn2_cl_ch3cl": (-1, 1, [("Cl", 0, 0, -3.20), ("C", 0, 0, 0), ("Cl", 0, 0, 1.80),
                             ("H", 1.03, 0, -0.36), ("H", -0.515, 0.892, -0.36), ("H", -0.515, -0.892, -0.36)],
                     [("Cl", 0, 0, -1.80), ("C", 0, 0, 0), ("Cl", 0, 0, 3.20),
                      ("H", 1.03, 0, 0.36), ("H", -0.515, 0.892, 0.36), ("H", -0.515, -0.892, 0.36)]),
    "h_hcl": (0, 2, [("H", 0, 0, -2.50), ("H", 0, 0, 0), ("Cl", 0, 0, 1.28)],
              [("H", 0, 0, -0.74), ("H", 0, 0, 0), ("Cl", 0, 0, 3.00)]),
    "oh_h2": (0, 2, [("O", 0, 0, 0), ("H", -0.94, 0.25, 0), ("H", 0, 0, 2.40), ("H", 0, 0, 3.14)],
              [("O", 0, 0, 0), ("H", -0.94, 0.25, 0), ("H", 0.24, 0, 0.94), ("H", 0, 0, 3.20)]),
    "ch3_h2": (0, 2, [("C", 0, 0, 0), ("H", 1.08, 0, 0), ("H", -0.54, 0.935, 0), ("H", -0.54, -0.935, 0),
                      ("H", 0, 0, 2.50), ("H", 0, 0, 3.24)],
               [("C", 0, 0, 0), ("H", 1.02, 0, -0.36), ("H", -0.51, 0.885, -0.36), ("H", -0.51, -0.885, -0.36),
                ("H", 0, 0, 1.09), ("H", 0, 0, 3.30)]),
}


def xyz(path, atoms, title):
    with open(path, "w") as fh:
        fh.write(f"{len(atoms)}\n{title}\n")
        for el, x, y, z in atoms:
            fh.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def write(name, base, method="NEB-TS", images=8, extra=""):
    charge, mult, rea, pro = R[base]
    d = os.path.join(OUT, name)
    os.makedirs(d, exist_ok=True)
    xyz(os.path.join(d, "reactant.xyz"), rea, f"{base} reactant")
    xyz(os.path.join(d, "product.xyz"), pro, f"{base} product")
    uks = "UKS " if mult > 1 else ""
    freq = " Freq" if method == "NEB-TS" else ""
    with open(os.path.join(d, f"{name}.inp"), "w") as fh:
        fh.write(f"! {uks}B3LYP D3BJ def2-SVP TightSCF {method}{freq}\n%neb\n  NEB_End_XYZFile \"product.xyz\"\n"
                 f"  NImages {images}\n  PreOpt_Ends true\n{extra}end\n* xyzfile {charge} {mult} reactant.xyz\n")


def main():
    for name in R:
        if name != "sn2_f_ch3cl":
            write(name, name)
    write("h2co_dissociation_3img", "h2co_dissociation", method="NEB-CI", images=3, extra="  Convtype CIonly\n")
    write("sn2_f_ch3cl_maxiter3", "sn2_f_ch3cl", extra="  MaxIter 3\n")
    write("oh_h2_maxiter3", "oh_h2", extra="  MaxIter 3\n")
    print(len(os.listdir(OUT)), "inputs")


if __name__ == "__main__":
    main()
