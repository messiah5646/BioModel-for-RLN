"""L1 measurement contracts (SYSTEM-DESIGN v2 §4, kept in v3).

An observable is not a float: it is a value plus its scale, dispersion type, per-visit n, population and source.
Each rule below is traceable to a near-miss in rounds 55-58 (see docs/SYSTEM-DESIGN-v3.md §4).
"""
from __future__ import annotations
from dataclasses import dataclass, field
import math

SCALES = {"absolute", "change_from_baseline", "ratio_to_baseline", "pct_change", "geometric_ratio", "responders"}
DISPERSIONS = {"sd", "se", "ci95", "iqr", "none"}
POPULATIONS = {"ITT", "mITT", "completers", "continuers", "responders", "per-protocol"}
KNOWN_MISSING = (99999.0, 9999999.0, 99.999)


class ContractError(ValueError):
    pass


@dataclass(frozen=True)
class Observable:
    value: float
    scale: str
    dispersion: str
    n: int                       # per visit AND per endpoint, never per trial
    population: str
    unit: str
    source: str                  # URL or file + retrieval date
    anchor: str | None = None    # which baseline a relative scale refers to; None => not convertible
    spread: float | None = None  # the dispersion number itself (SD, SE, CI half-width...)
    missing_code: float | None = None

    def __post_init__(self):
        if self.scale not in SCALES:
            raise ContractError(f"unknown scale {self.scale!r}")
        if self.dispersion not in DISPERSIONS:
            raise ContractError(f"unknown dispersion {self.dispersion!r}")
        if self.population not in POPULATIONS:
            raise ContractError(f"unknown population {self.population!r}")
        if self.n <= 0:
            raise ContractError("n must be the per-visit, per-endpoint analysed count (> 0)")
        if self.value in KNOWN_MISSING or (self.missing_code is not None and self.value == self.missing_code):
            raise ContractError(f"value {self.value} is a missing-data code (rule: missing codes are never numbers)")
        if self.scale in {"change_from_baseline", "ratio_to_baseline", "pct_change"} and self.anchor is None:
            # allowed to exist, but marked unconvertible (CFZ533 case, round 58)
            pass

    def sd(self) -> float:
        """Explicit, recorded dispersion conversion (rule 1). Never implicit."""
        if self.spread is None or self.dispersion == "none":
            raise ContractError("no dispersion recorded")
        if self.dispersion == "sd":
            return self.spread
        if self.dispersion == "se":
            return self.spread * math.sqrt(self.n)
        if self.dispersion == "ci95":
            return self.spread / 1.96 * math.sqrt(self.n)
        raise ContractError(f"{self.dispersion} cannot be converted to an SD without distributional assumptions")


def combinable(a: Observable, b: Observable) -> None:
    """Raise unless two observables may be differenced, pooled or concatenated (rules 2-4)."""
    if a.scale != b.scale:
        raise ContractError(f"cross-scale arithmetic: {a.scale} vs {b.scale}")
    if a.population != b.population:
        raise ContractError(f"cross-population concatenation: {a.population} vs {b.population}")
    if a.unit != b.unit:
        raise ContractError(f"unit mismatch: {a.unit} vs {b.unit} (unit fields are checked, not trusted)")
    if a.scale != "absolute" and a.anchor != b.anchor:
        raise ContractError("relative scales with different anchors are not comparable")


@dataclass(frozen=True)
class ArmCount:
    """Responders in a drug arm vs its control arm at one visit: the unit the drug-effect layer consumes."""
    trial: str
    arm: str
    week: float
    r_drug: int
    n_drug: int
    r_ctrl: int
    n_ctrl: int
    population: str = "ITT"
    source: str = ""

    def __post_init__(self):
        for r, n in ((self.r_drug, self.n_drug), (self.r_ctrl, self.n_ctrl)):
            if not (0 <= r <= n) or n <= 0:
                raise ContractError(f"{self.trial}/{self.arm}: responders {r} outside 0..{n}")

    def log_or(self, cc: float = 0.5) -> tuple[float, float]:
        """Log odds ratio and SE with the project's 0.5 continuity correction (same as lumen_forecast.load)."""
        a, b = self.r_drug + cc, self.n_drug - self.r_drug + cc
        c, d = self.r_ctrl + cc, self.n_ctrl - self.r_ctrl + cc
        return math.log(a * d / (b * c)), math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)


def check_visits(rows: list[ArmCount]) -> list[str]:
    """Rule 5: per-arm n must not rise across visits (cohort change) -> raise.
    Rule 6: identical counts at two visits -> returned as a warning. For continuous means (NOBILITY, round 58) this is a
    duplicate; for small integer counts it can be coincidence, so the source must be checked by hand, not auto-dropped."""
    warnings: list[str] = []
    by_arm: dict[tuple, list[ArmCount]] = {}
    for r in rows:
        by_arm.setdefault((r.trial, r.arm), []).append(r)
    for key, rs in by_arm.items():
        rs = sorted(rs, key=lambda x: x.week)
        for p, q in zip(rs, rs[1:]):
            if q.n_drug > p.n_drug or q.n_ctrl > p.n_ctrl:
                raise ContractError(f"{key}: n rises from week {p.week} to {q.week} (cohort change, not attrition)")
            if (p.r_drug, p.n_drug, p.r_ctrl, p.n_ctrl) == (q.r_drug, q.n_drug, q.r_ctrl, q.n_ctrl):
                warnings.append(f"{key}: weeks {p.week} and {q.week} have identical counts; check the source")
    return warnings
