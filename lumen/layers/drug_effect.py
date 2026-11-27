"""Drug-effect layer: relative effect (log OR, drug vs control) of a lupus/LN add-on drug.

The maths is experiments/lumen-bench/lumen_forecast.py, loaded from its md5-frozen file. It is imported, never copied,
so frozen forecasts (r90, W1 LN, W1 SLE) and audits keep reproducing. This module adds only the v3 output policy:
  * no --early: headline = pooled prior (round 116 / 127; Supported, n 50 trials time-forward)
  * with early read-out: NO headline; null, pooled and carry-forward side by side (rule 2b, round 144)
  * control rate: pooled, not modelled; no absolute arm rates (round 129)
"""
from __future__ import annotations
import importlib.util, math
from ..provenance import verify

TABLES = {
    "ln": "experiments/lumen-bench/resp_table.csv",          # LN complete renal response, 20 trials
    "ln_plus_pub": "experiments/lumen-bench/resp_table_plus_pub.csv",
    "sri4": "experiments/lumen-bench/sle/sri4_table.csv",    # non-renal SLE, SRI-4
    "bicla": "experiments/lumen-bench/sle/bicla_table.csv",  # non-renal SLE, BICLA
}
_CORE = None


def core():
    """The frozen forecaster module (lumen_forecast.py), md5-checked on first load."""
    global _CORE
    if _CORE is None:
        p = verify("experiments/lumen-bench/lumen_forecast.py")
        spec = importlib.util.spec_from_file_location("lumen_forecast_frozen", p)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _CORE = mod
    return _CORE


def table(name: str):
    return core().load(str(verify(TABLES[name])))


def prior(molecule: str, target: str, table_name: str = "ln", holdout: str | None = None, anchor: bool = False) -> dict:
    """Prior-only forecast. anchor=False is the shipped default (pooled); anchor=True reproduces frozen M3 / r90 maths."""
    return core().forecast(molecule, target, None, R=table(table_name), holdout=holdout, anchor=anchor)


def forecast(molecule: str, target: str, table_name: str = "ln", early=None, holdout: str | None = None,
             anchor: bool = False, use_early_update: bool = False) -> dict:
    """Same JSON as `python lumen_forecast.py MOL TARGET [--early ...]`, plus a v3 'layer' block."""
    L = core()
    R = table(table_name)
    out = L.forecast(molecule, target, None, R=R, holdout=holdout, anchor=anchor)
    if early:
        upd = L.forecast(molecule, target, tuple(early), R=R, holdout=holdout, anchor=anchor)
        if use_early_update:
            out = upd
        else:
            out["early_update_EMERGING"] = {k: upd[k] for k in ("early_logOR", "post_logOR", "post_sd", "OR", "OR_95") if k in upd}
            wk, rd, nd, rc, nc = early
            a, b, c, d = rd + .5, nd - rd + .5, rc + .5, nc - rc + .5
            le = math.log(a * d / (b * c)); se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
            out["early_note"] = ("No early read-out rule is validated: RT-1 (round 127) shrunk update lost to carry-forward by 77% "
                                 "(Emerging, n 7 trials); RT-3 and RT-4 UNDERPOWERED (rounds 136, 141); settled under rule 8 until "
                                 "W1 prospective results. Full text: python experiments/lumen-bench/lumen_forecast.py ... --early")
            out["headline"] = ("NONE: no validated rule; pooled lost to carry-forward on the all-trial early subsets "
                               "(RT-1, RT-4) and to null in LN")
            out["early_side_by_side"] = {
                "M0_null": {"OR": 1.0},
                "M1_pooled_prior": {"OR": out["OR"], "OR_95": out["OR_95"]},
                "C0_carry_forward": {"OR": round(math.exp(le), 2),
                                     "OR_95_sampling_only": [round(math.exp(le - 1.96 * se), 2), round(math.exp(le + 1.96 * se), 2)],
                                     "note": "raw early read-out (week %g) carried to the late visit; early->late drift not added" % wk}}
    out["control_rate"] = "pooled, not modelled; absolute arm rates not output (round 129, RT-2 PREREG d5b5f84)"
    out["layer"] = {"name": "drug_effect", "table": table_name, "status": "shipped",
                    "evidence": "pooled prior beats null time-forward: Supported, n 50 trials (round 127)"}
    return out
