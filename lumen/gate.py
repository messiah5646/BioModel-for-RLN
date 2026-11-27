"""L2 gate: the rule every module must pass before it ships (queue rules 2, 2b, 12).

verdict(candidate, baselines, ...) returns one of
  UNDERPOWERED  fewer held-out units than the PREREG minimum
  FAIL          loses to ANY trivial baseline by >= margin (rule 2b: a real FAIL branch, not 'inconclusive')
  PASS          beats the reference baseline by >= margin and loses to none
  INCONCLUSIVE  everything else
Lower scores are better (IV-MSE, log score). This is the same decision logic as timeforward/run_rt1.summarise_early,
generalised so that new modules (patient-state, cell-state, a tau prior) are judged the same way.
"""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class GateResult:
    verdict: str
    candidate: float
    baselines: dict
    reference: str
    n_units: int
    n_trials: int
    margin: float
    reasons: list = field(default_factory=list)

    def ships(self) -> bool:
        return self.verdict == "PASS"


def verdict(candidate: float, baselines: dict, reference: str, n_units: int, n_trials: int,
            min_units: int = 6, min_trials: int = 4, margin: float = 0.10) -> GateResult:
    if reference not in baselines:
        raise ValueError("the reference baseline must be one of the baselines")
    r = GateResult("INCONCLUSIVE", candidate, dict(baselines), reference, n_units, n_trials, margin)
    if n_units < min_units or n_trials < min_trials:
        r.verdict = "UNDERPOWERED"
        r.reasons.append(f"n {n_units} units / {n_trials} trials < minimum {min_units} / {min_trials}")
        return r
    best_name = min(baselines, key=baselines.get)
    if candidate >= (1 + margin) * baselines[best_name]:
        r.verdict = "FAIL"
        r.reasons.append(f"loses to {best_name} by {candidate / baselines[best_name] - 1:.0%} (rule 2b)")
        return r
    beats_ref = candidate <= (1 - margin) * baselines[reference]
    loses_none = all(candidate < (1 + margin) * v for v in baselines.values())
    if beats_ref and loses_none:
        r.verdict = "PASS"
        r.reasons.append(f"beats {reference} by {1 - candidate / baselines[reference]:.0%} and loses to no baseline")
    else:
        r.reasons.append("does not clear the margin over the reference")
    return r
