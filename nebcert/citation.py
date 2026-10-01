"""
Citation strings for NEBCert (concept DOI: resolves to the latest version on Zenodo).
"""

from nebcert import __version__

CONCEPT_DOI = "10.5281/zenodo.22217586"
TITLE = "NEBCert: quality checks for NEB paths, transition states and tunnelling-corrected rate constants"

BIBTEX = (
    "@software{monreal2026nebcert,\n"
    "  author = {Monreal-Hern{\\'a}ndez, Andr{\\'e}s},\n"
    f"  title = {{{{{TITLE}}}}},\n"
    "  year = {2026},\n"
    f"  version = {{{__version__}}},\n"
    "  publisher = {Zenodo},\n"
    f"  doi = {{{CONCEPT_DOI}}},\n"
    "  url = {https://github.com/sircalch/nebcert}\n"
    "}\n"
)

APA = f"Monreal-Hernández, A. (2026). {TITLE} (v{__version__}). Zenodo. https://doi.org/{CONCEPT_DOI}"
