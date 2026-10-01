"""
Aggregation of the NEB, transition-state, rate-constant and tunnelling checks into one report.

The overall status is the worst individual status (FAIL > WARNING > PASS); a rate constant that the input
does not determine (electronic barrier) is reported as NOT_APPLICABLE and does not count. A PASS means that the
checks run found nothing wrong; it does not establish that the path or the barrier is correct.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
import numpy as np

from nebcert import __version__

from nebcert.core.neb_profile import NEBProfileResult
from nebcert.core.ts_frequency import TSFrequencyResult
from nebcert.core.tst_kinetics import TSTKineticsResult
from nebcert.core.tunneling import TunnelingResult


@dataclass
class ReactionPathwayReport:
    overall_status: str  # 'PASS', 'WARNING', 'FAIL', 'NOT_APPLICABLE'
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
        # NOT_APPLICABLE (rate from an electronic barrier) is reported but does not set the overall status
        if tst_res.status != "NOT_APPLICABLE":
            statuses.append(tst_res.status)
        if tst_res.status != "PASS":
            recommendations.append(tst_res.diagnostic_message)

    if tunneling_res is not None:
        statuses.append(tunneling_res.status)
        if tunneling_res.status != "PASS":
            recommendations.append(tunneling_res.diagnostic_message)

    if not statuses:
        overall_status = "NOT_APPLICABLE"
        validation_score = "NO CHECKS RUN"
    elif "FAIL" in statuses:
        overall_status = "FAIL"
        validation_score = "AT LEAST ONE CHECK FAILED"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        validation_score = "PASSED WITH WARNINGS"
    else:
        overall_status = "PASS"
        validation_score = "ALL CHECKS PASSED"

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
            "version": __version__,
            "citation": f"Monreal-Hernández, A. (2026). NEBCert: quality checks for NEB paths, transition states and tunnelling-corrected rate constants (v{__version__})."
        }
    )
