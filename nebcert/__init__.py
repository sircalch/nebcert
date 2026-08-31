"""
NEBCert: Automated Quality-Control, Transition State Verification,
Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification.
"""

__version__ = "1.0.0"
__author__ = "Andre Monreal-Hernández"
__license__ = "MIT"

from nebcert.core.neb_profile import calculate_neb_profile_analysis, NEBProfileResult
from nebcert.core.ts_frequency import verify_ts_frequency_and_irc, TSFrequencyResult
from nebcert.core.tst_kinetics import calculate_eyring_tst_rates, TSTKineticsResult
from nebcert.core.tunneling import calculate_quantum_tunneling_corrections, TunnelingResult
from nebcert.core.scoring import assess_reaction_pathway_quality, ReactionPathwayReport

__all__ = [
    "__version__",
    "calculate_neb_profile_analysis",
    "NEBProfileResult",
    "verify_ts_frequency_and_irc",
    "TSFrequencyResult",
    "calculate_eyring_tst_rates",
    "TSTKineticsResult",
    "calculate_quantum_tunneling_corrections",
    "TunnelingResult",
    "assess_reaction_pathway_quality",
    "ReactionPathwayReport"
]
