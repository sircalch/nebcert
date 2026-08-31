"""
Multi-metric kinetics certification scoring, reaction pathway assessment, and report aggregation for NEBCert.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
import numpy as np

from nebcert.core.neb_profile import NEBProfileResult
from nebcert.core.ts_frequency import TSFrequencyResult
from nebcert.core.tst_kinetics import TSTKineticsResult
from nebcert.core.tunneling import TunnelingResult


@dataclass
class ReactionPathwayReport:
    overall_status: str  # 'PASS', 'WARNING', 'FAIL'
    validation_score: str
    metadata: Dict[str, Any]
    neb_profile: Optional[NEBProfileResult]
    ts_frequency: Optional[TSFrequencyResult]
    tst_kinetics: Optional[TSTKineticsResult]
    tunneling: Optional[TunnelingResult]
    recommendations: List[str]
    provenance: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_reaction_pathway_quality(
    metadata: Dict[str, Any],
    neb_res: Optional[NEBProfileResult] = None,
    ts_freq_res: Optional[TSFrequencyResult] = None,
    tst_res: Optional[TSTKineticsResult] = None,
    tunneling_res: Optional[TunnelingResult] = None
) -> ReactionPathwayReport:
    """
    Consolidates NEB minimum energy path, TS frequency verification, TST kinetics,
    and quantum tunneling corrections.

    Parameters
    ----------
    metadata : dict
        Reaction title, software engine, exchange-correlation functional / level of theory.
    neb_res : NEBProfileResult, optional
    ts_freq_res : TSFrequencyResult, optional
    tst_res : TSTKineticsResult, optional
    tunneling_res : TunnelingResult, optional

    Returns
    -------
    report : ReactionPathwayReport
    """
    statuses = []
    recommendations = []

    if neb_res is not None:
        statuses.append(neb_res.status)
        if neb_res.status != "PASS":
            recommendations.append(neb_res.diagnostic_message)

    if ts_freq_res is not None:
        statuses.append(ts_freq_res.status)
        if ts_freq_res.status != "PASS":
            recommendations.append(ts_freq_res.diagnostic_message)

    if tst_res is not None:
        statuses.append(tst_res.status)
        if tst_res.status != "PASS":
            recommendations.append(tst_res.diagnostic_message)

    if tunneling_res is not None:
        statuses.append(tunneling_res.status)
        if tunneling_res.status != "PASS":
            recommendations.append(tunneling_res.diagnostic_message)

    if not statuses:
        overall_status = "PASS"
        validation_score = "REACTION PATHWAY AUDIT = UNVERIFIED"
    elif "FAIL" in statuses:
        overall_status = "FAIL"
        validation_score = "REACTION PATHWAY AUDIT = FAILED / UNPHYSICAL ARTIFACTS DETECTED"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        validation_score = "REACTION PATHWAY AUDIT = ACCEPTABLE WITH METHODOLOGICAL WARNINGS"
    else:
        overall_status = "PASS"
        validation_score = "REACTION PATHWAY AUDIT = FULLY VALIDATED (PUBLICATION GRADE)"

    return ReactionPathwayReport(
        overall_status=overall_status,
        validation_score=validation_score,
        metadata=metadata,
        neb_profile=neb_res,
        ts_frequency=ts_freq_res,
        tst_kinetics=tst_res,
        tunneling=tunneling_res,
        recommendations=recommendations,
        provenance={
            "tool": "NEBCert",
            "version": "1.0.0",
            "citation": "Monreal-Hernández, A. (2026). NEBCert: Automated Quality-Control, Transition State Verification, Nudged Elastic Band (NEB), Quantum Tunneling, and Reaction Kinetics Certification."
        }
    )
