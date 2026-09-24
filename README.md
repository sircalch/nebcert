# NEBCert

[![CI](https://github.com/sircalch/nebcert/actions/workflows/test.yml/badge.svg)](https://github.com/sircalch/nebcert/actions)
[![PyPI version](https://img.shields.io/pypi/v/nebcert.svg?color=blue)](https://pypi.org/project/nebcert/)
[![Python versions](https://img.shields.io/pypi/pyversions/nebcert.svg)](https://pypi.org/project/nebcert/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234600.svg)](https://doi.org/10.5281/zenodo.1234600)

> **Automated Quality-Control, Transition State Verification, Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification (VASP, ORCA, Gaussian).**

---

## Overview

**NEBCert** is an open-source scientific software package designed to standardize, audit, and certify chemical reaction pathways, Nudged Elastic Band (NEB / CI-NEB) minimum energy paths (MEPs), transition state 1st-order saddle points, Eyring-Polanyi TST rate constants $k(T)$, and quantum mechanical tunneling corrections (Wigner & asymmetric Eckart barrier transmission).

In computational chemistry, heterogeneous catalysis, and atmospheric/combustion kinetics:

- 🧗 **Minimum Energy Path (MEP) & CI-NEB Spline Analysis**:
  - Fits natural cubic splines along normalized cumulative reaction coordinates $s$.
  - Extracts accurate forward barrier $E_a^{\text{fwd}}$, reverse barrier $E_a^{\text{rev}}$, and reaction energy $\Delta E_{\text{rxn}}$ (in eV and kcal/mol).
  - Audits climbing image tangent residual forces ($|F_\parallel| \le 0.05\text{ eV/\AA}$) and detects unphysical intermediate kinks / hidden intermediates.
- 🎯 **Transition State (TS) & IRC Audit**:
  - Certifies true 1st-order saddle points (strictly 1 imaginary vibrational frequency $\nu_{\text{imag}} < 0\text{ cm}^{-1}$ along the reaction bond coordinate).
  - Detects artifactual 2nd-order saddle points ($> 1$ imaginary frequencies) or false local minima (0 imaginary modes).
  - Verifies Intrinsic Reaction Coordinate (IRC) connectivity connecting reactants to products.
- ⏱️ **Eyring Transition State Theory (TST) & Arrhenius Kinetics**:
  - $k_{\text{TST}}(T) = \frac{k_B T}{h} \exp\left(-\frac{\Delta G^\ddagger(T)}{R T}\right)$ across user-defined temperature ranges (200 to 1000 K).
  - Fits pre-exponential factor $A$ ($\text{s}^{-1}$) and Arrhenius activation energy $E_a$ ($\text{kcal/mol}$).
- 🌌 **Quantum Mechanical Tunneling Corrections**:
  - **Wigner Tunneling Factor**: $\kappa_W(T) = 1 + \frac{1}{24} \left(\frac{h |\nu^\ddagger|}{k_B T}\right)^2$.
  - **Asymmetric Eckart Barrier Transmission**: Numerical Gauss-Legendre quadrature integration of 1D asymmetric Eckart potential $P(E)$, predicting low-temperature quantum tunneling accelerations.
- 📑 **Publication Deliverables**:
  - Interactive self-contained `report.html` dashboard.
  - Publication vector figures (NEB MEP Spline Profile, Arrhenius Plot $\ln k$ vs $1/T$, Tunneling $\kappa(T)$) in SVG, PDF, PNG (300 DPI).
  - Ready-to-compile LaTeX summary tables (`.tex`).
  - Draft **Methods** text snippet and BibTeX citation (`citation.bib`).

```
       NEB Outputs (VASP neb.dat, ORCA NEB, Gaussian IRC, CSV)
                               │
                               ▼
  ┌───────────────────────────────────────────────────────────┐
  │                          NEBCert                          │
  │  ├── NEB / CI-NEB Minimum Energy Path Spline Interpolation │
  │  ├── TS 1st-Order Saddle Point & Imaginary Freq Verification│
  │  ├── Eyring Transition State Theory Rate Constants k(T)   │
  │  └── Asymmetric Eckart & Wigner Quantum Tunneling Factors │
  └───────────────────────────────────────────────────────────┘
                               │
                               ▼
  ┌───────────────────────────────────────────────────────────┐
  │                   Publication Deliverables                │
  │  ├── report.html (Interactive Dashboard & Badges)         │
  │  ├── nebcert_neb_mep_profile.pdf/svg/png                  │
  │  ├── nebcert_arrhenius_kinetics.pdf/svg/png               │
  │  ├── nebcert_tunneling_factors.pdf/svg/png                │
  │  ├── nebcert_summary_table.tex / .csv                     │
  │  ├── methods_snippet.txt (Ready for Manuscript)           │
  │  └── citation.bib (BibTeX Reference)                      │
  └───────────────────────────────────────────────────────────┘
```

---

## Installation

### From PyPI
> **Note:** PyPI release pending. Until then, install from the tagged GitHub release:

```bash
pip install "git+https://github.com/sircalch/nebcert@v1.0.0"
```

### From Source
```bash
git clone https://github.com/sircalch/nebcert.git
cd nebcert
pip install -e .[dev]
```

---

## Quickstart (CLI)

### 1. Run Benchmark Demo (H-abstraction 9-Image NEB + TS + TST + Tunneling)
```bash
nebcert demo -o my_kinetics_audit/
```
Open `my_kinetics_audit/report.html` in any browser!

### 2. Assess NEB neb.dat File
```bash
nebcert assess -i neb.dat --frequencies "-1250,150,300,800" --reaction "CH4 + OH -> CH3 + H2O" -o neb_report/
```

---

## Python API Usage

```python
from nebcert import (
    calculate_neb_profile_analysis,
    verify_ts_frequency_and_irc,
    calculate_eyring_tst_rates,
    calculate_quantum_tunneling_corrections,
    assess_reaction_pathway_quality
)
from nebcert.reporters import (
    generate_nebcert_figures,
    generate_nebcert_manuscript_assets,
    generate_nebcert_html_report
)

# 1. NEB Profile Analysis
neb_res = calculate_neb_profile_analysis(
    energies_ev=[0.0, 0.05, 0.15, 0.23, 0.10, -0.25, -0.62],
    coordinates_s_ang=[0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]
)

# 2. Transition State Verification
ts_res = verify_ts_frequency_and_irc(
    frequencies_cm1=[-1250.0, 120.0, 300.0, 1500.0],
    irc_confirmed=True
)

# 3. Eyring TST Kinetics & Quantum Tunneling
tst_res = calculate_eyring_tst_rates(e_activation_ev=neb_res.e_forward_barrier_ev)
tun_res = calculate_quantum_tunneling_corrections(
    imaginary_freq_cm1=1250.0,
    e_forward_barrier_ev=neb_res.e_forward_barrier_ev,
    e_reverse_barrier_ev=neb_res.e_reverse_barrier_ev
)

# 4. Consolidate reaction pathway report
report = assess_reaction_pathway_quality(
    metadata={"reaction": "CH4 + OH -> CH3 + H2O", "functional": "wB97X-D3", "software": "ORCA"},
    neb_res=neb_res,
    ts_freq_res=ts_res,
    tst_res=tst_res,
    tunneling_res=tun_res
)

print(f"Overall Certification: {report.overall_status}")
print(f"Forward Barrier: {report.neb_profile.e_forward_barrier_kcal_mol:.2f} kcal/mol")
print(f"TST Rate (298 K): {report.tst_kinetics.k_298_s_minus_1:.2e} s^-1")
print(f"Tunneling Factor: kappa_Eckart = {report.tunneling.kappa_eckart_298:.2f}")

# 5. Export deliverables
generate_nebcert_figures(report, "output_dir/")
generate_nebcert_manuscript_assets(report, "output_dir/")
generate_nebcert_html_report(report, "output_dir/report.html")
```

---

## Citation

If you use NEBCert in your publications, please cite:

```bibtex
@software{monreal2026nebcert,
  author = {Monreal-Hern{\'a}ndez, Andre},
  title = {{NEBCert: Automated Quality-Control, Transition State Verification, Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/nebcert}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

