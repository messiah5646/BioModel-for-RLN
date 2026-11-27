"""L4 module registry: what LUMEN ships, at which evidence label and n, and what it would take to change that.

A module may produce headline output only if status == 'shipped', which requires a gate verdict of PASS on a
pre-registered held-out test (lumen.gate). Everything else is printed, if at all, as 'not validated'.
Entries copy claude/where-this-stands.md (round 145); they add no new claims.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict

STATUSES = ("shipped", "option-emerging", "standalone", "failed", "not-implemented", "supporting")


@dataclass(frozen=True)
class Module:
    name: str
    layer: str
    status: str
    evidence: str          # label + n, as in where-this-stands
    rounds: str
    consumes: str
    produces: str
    reopen_if: str         # the one thing that would change the status

    def __post_init__(self):
        if self.status not in STATUSES:
            raise ValueError(self.status)
        if self.status == "shipped" and "n " not in self.evidence:
            raise ValueError(f"{self.name}: a shipped module must state n (rule 12)")


REGISTRY = [
    Module("drug_effect.pooled_prior", "drug-effect", "shipped",
           "Supported, n 50 trials time-forward (LN 15, SLE 35); beats null; anchoring adds nothing (n 23 SLE trials)",
           "116, 127, 141", "arm counts (resp/sri4/bicla tables)", "log OR + 95% interval, P(benefit)",
           "W1 prospective scoring (>= 5 frozen SLE trials) under PREREG-W1-SLE"),
    Module("drug_effect.anchor_molecule_target", "drug-effect", "failed",
           "anchor 31% worse than pooled, Supported n 23 trials (SLE); LN missed margin (20 trials)",
           "86, 116", "arm counts", "log OR", "kept only as frozen M3 secondary"),
    Module("drug_effect.early_update", "drug-effect", "option-emerging",
           "lost to carry-forward by 77%, Emerging n 7 trials (RT-1); RT-3/RT-4 UNDERPOWERED; settled (rule 8)",
           "112, 127, 136, 141", "early by-arm read-out", "posterior log OR (flag only)",
           "W1 prospective results, >= 5 frozen trials reporting both visits"),
    Module("control_rate", "drug-effect", "not-implemented",
           "meta-regression did not beat pooled: INCONCLUSIVE 1.06x, n 60 trials; LN FAIL 1.17x, n 20",
           "129", "control arms", "none (pooled, not modelled)", "a new covariate pre-registered against the pooled rate"),
    Module("patient_state.ifn_setpoint", "patient-state", "standalone",
           "IFN set point Supported n 135 placebo patients; only 6-11% better than carry-forward; no consumer",
           "96, 117, 121, 123", "baseline IFN score", "IFN(t)", "a LUMEN module that consumes f(t) before day 90"),
    Module("cell_state.vc1_learned", "cell-state", "failed",
           "LUMEN-JEPA (PCA+MLP) 1.32x worse than population average, 0/11 wins, Emerging n 11",
           "142, 143", "baseline pseudobulk", "month-3 pseudobulk",
           "an open longitudinal sc set with >= 30 patients (W2 watches)"),
    Module("flare", "flare", "not-implemented",
           "baseline BCMA/IFN/SLEDAI do not order time to flare, Supported negative n 39 patients",
           "122", "clinical history", "none (clinical-only)", "new data type"),
    Module("mechanistic_core", "mechanistic", "supporting",
           "loses to 0-1 parameter rules on 6 observables (rounds 50-85); used for exchangeability and horizons only",
           "50, 53, 61, 78, 82, 83, 85", "parameters", "identifiability verdicts", "beats persistence out-of-sample on levels"),
]


def shipped():
    return [m for m in REGISTRY if m.status == "shipped"]


def as_table():
    return [asdict(m) for m in REGISTRY]
