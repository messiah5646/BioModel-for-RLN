"""REC1 freeze (PREREG-REC1.md, 6f95104). Write-once prospective freeze of three recalibration competitors
(a) mixture-weight, (b) van Zwet SNR shrinkage, (c) equal-weight median ensemble, against B (printed base rate),
for the 21 W1 (LN) + W1-SLE trials. No competitor is scored on historical trials (rule 8): cal1_units.csv is used
only to fit w (and the descriptive temperature T).

Usage:  python experiments/lumen-v3/freeze_rec1.py --selftest     (simulated inputs only; writes nothing in the repo)
        python experiments/lumen-v3/freeze_rec1.py                (real freeze; write-once)
"""
from __future__ import annotations
import argparse, hashlib, json, math, pathlib, sys
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
PREREG = "experiments/lumen-v3/PREREG-REC1.md"
CAL1_UNITS = REPO / "experiments/lumen-v3/cal1_units.csv"
CAL1_UNITS_MD5 = "6243f7969c75b7fd66e6bfee39a8ec3d"          # = git show 2b11bca:experiments/lumen-v3/cal1_units.csv
PROSP = REPO / "experiments/lumen-bench/prospective"
W1_LN = PROSP / "w1/w1_frozen_2026-11-21.json"
W1_SLE = PROSP / "w1-sle/frozen_2026-11-26.json"
SRC = PROSP / "rec1/sources"
OUT_DIR = PROSP / "rec1"

# van Zwet, Schwab & Senn 2021 (Stat Med, doi 10.1002/sim.9173), Table 1: zero-mean normal mixture for the SNR
VZ_W = np.array([0.32, 0.31, 0.30, 0.07])
VZ_TAU = np.array([0.61, 1.42, 2.16, 5.64])
CONTROL_TYPES = {"PLACEBO_COMPARATOR", "SHAM_COMPARATOR", "SHAM", "NO_INTERVENTION"}
TABLE_KEY = {"sle/sri4_table.csv": "sri4", "sle/bicla_table.csv": "bicla"}


# ---------------------------------------------------------------- competitors (pure functions, tested in selftest)
def se_ant2(n_a: float, p0: float, M1: float) -> float:
    p1 = 1 / (1 + math.exp(-(math.log(p0 / (1 - p0)) + M1)))
    return 1 / (n_a * p1) + 1 / (n_a * (1 - p1)) + 1 / (n_a * p0) + 1 / (n_a * (1 - p0))


def P1_of(M1, v, s2):
    return float(norm.cdf(M1 / math.sqrt(v + s2)))


def P_a(P1, w):
    return (1 - w) * P1 + w * 0.5


def P_T(P1, T):
    return float(norm.cdf(norm.ppf(P1) / T))


def vz_posterior(z):
    """Posterior of SNR | z under the van Zwet mixture: weights q_k, means m_k, variances s2_k."""
    t2 = VZ_TAU ** 2
    lik = VZ_W * norm.pdf(z, 0, np.sqrt(t2 + 1))
    q = lik / lik.sum()
    return q, z * t2 / (t2 + 1), t2 / (t2 + 1)


def P_b(M1, v, s2):
    sd = math.sqrt(v)
    q, m, s2k = vz_posterior(M1 / sd)
    return float(np.sum(q * norm.cdf(sd * m / np.sqrt(v * s2k + s2)))), (q, m, s2k)


def mixture_interval(weights, means, sds, level):
    """Central interval of a normal mixture (numerical quantiles)."""
    weights, means, sds = map(np.asarray, (weights, means, sds))
    F = lambda x: float(np.sum(weights * norm.cdf((x - means) / sds)))
    lo_p, hi_p = (1 - level) / 2, 1 - (1 - level) / 2
    a, b = float(np.min(means - 12 * sds)), float(np.max(means + 12 * sds))
    return [round(brentq(lambda x: F(x) - lo_p, a, b), 4), round(brentq(lambda x: F(x) - hi_p, a, b), 4)]


def fit_w(P1, E, wt):
    grid = np.round(np.arange(0, 1.0005, 0.001), 3)
    ll = [np.sum(wt * (E * np.log(P_a(P1, w)) + (1 - E) * np.log(1 - P_a(P1, w)))) for w in grid]
    return float(grid[int(np.argmax(ll))])


def fit_T(P1, E, wt):
    grid = np.round(np.arange(0.50, 5.0005, 0.01), 2)
    z = norm.ppf(np.clip(P1, 1e-12, 1 - 1e-12))
    def ll(T):
        p = np.clip(norm.cdf(z / T), 1e-12, 1 - 1e-12)
        return np.sum(wt * (E * np.log(p) + (1 - E) * np.log(1 - p)))
    return float(grid[int(np.argmax([ll(T) for T in grid]))])


def trial_forecast(M1, v, s2, B, w, T, cf=None):
    """All frozen numbers for one trial. cf = (M3, v3) if the anchor is the same molecule, else None."""
    P1 = P1_of(M1, v, s2)
    pa = P_a(P1, w)
    pb, (q, m, s2k) = P_b(M1, v, s2)
    P_cf = P1_of(cf[0], cf[1], s2) if cf is not None else B
    pc = float(np.median([P1, 0.5, P_cf]))
    sdp = math.sqrt(v + s2); sd = math.sqrt(v)
    iv = {}
    for lev in (0.5, 0.95):
        k = "%d" % round(lev * 100)
        iv["P1_and_c_member_" + k] = [round(M1 - norm.ppf(.5 + lev / 2) * sdp, 4), round(M1 + norm.ppf(.5 + lev / 2) * sdp, 4)]
        iv["a_" + k] = mixture_interval([1 - w, w], [M1, 0.0], [sdp, math.sqrt(s2)], lev)
        iv["b_" + k] = mixture_interval(q, sd * m, np.sqrt(v * s2k + s2), lev)
    return {"P1": round(P1, 4), "P_a": round(pa, 4), "P_b": round(pb, 4), "P_c": round(pc, 4), "P_cf_member": round(P_cf, 4),
            "B": B, "D_legacy_latent": round(float(norm.cdf(M1 / sd)), 4), "C05": 0.5,
            "P_T_descriptive": round(P_T(P1, T), 4),
            "vz_posterior": {"z": round(M1 / sd, 4), "q": [round(x, 4) for x in q], "m": [round(x, 4) for x in m]},
            "intervals_logOR": iv}


# ---------------------------------------------------------------- inputs
def arms_from_ctgov(rec):
    arms = rec["protocolSection"].get("armsInterventionsModule", {}).get("armGroups", [])
    ctrl = [a["label"] for a in arms if a.get("type") in CONTROL_TYPES
            or "placebo" in a["label"].lower() or "standard of care" in a["label"].lower()]
    return len(arms), ctrl


def table_p0(R):
    P = R[R.week <= 104].sort_values("week").groupby(["trial", "arm"]).tail(1)
    per_trial = (P.r_ctrl / P.n_ctrl).groupby(P.trial).mean()
    return float(np.median(per_trial)), int(len(per_trial))


def md5(p):
    return hashlib.md5(pathlib.Path(p).read_bytes()).hexdigest()


def fit_on_cal1():
    import pandas as pd
    assert md5(CAL1_UNITS) == CAL1_UNITS_MD5, "cal1_units.csv drifted from 2b11bca"
    U = pd.read_csv(CAL1_UNITS)
    wt = 1 / U.groupby("trial").arm.transform("count").values
    P1, E = U.P1.values, U.E1.values
    return fit_w(P1, E, wt), fit_T(P1, E, wt), int(len(U)), int(U.trial.nunique())


# ---------------------------------------------------------------- selftest (simulated inputs only)
def selftest():
    rng = np.random.default_rng(1)
    # 1. w recovery on simulated arms
    P1 = rng.uniform(0.55, 0.99, 4000); w_true = 0.4
    E = (rng.uniform(size=P1.size) < P_a(P1, w_true)).astype(float)
    w_hat = fit_w(P1, E, np.ones_like(P1)); assert abs(w_hat - w_true) < 0.08, w_hat
    T_hat = fit_T(P1, E, np.ones_like(P1)); assert T_hat > 1, T_hat
    # 2. van Zwet: z = 0 -> 0.5; shrinks toward 0.5 vs P1 for moderate z; posterior weights sum 1
    assert abs(P_b(0.0, 0.09, 0.2)[0] - 0.5) < 1e-12
    q, m, s2k = vz_posterior(1.5); assert abs(q.sum() - 1) < 1e-12 and np.all(np.abs(m) < 1.5)
    pb = P_b(0.3, 0.04, 0.1)[0]; assert 0.5 < pb < P1_of(0.3, 0.04, 0.1), pb
    # 3. mixture interval of a single normal equals the analytic interval
    lo, hi = mixture_interval([1.0], [0.2], [0.5], 0.95); assert abs(lo - (0.2 - 1.96 * 0.5)) < 1e-3 and abs(hi - (0.2 + 1.96 * 0.5)) < 1e-3
    # 4. se_ant: symmetric, decreasing in n
    assert se_ant2(100, 0.3, 0.3) < se_ant2(50, 0.3, 0.3)
    # 5. median ensemble with no same-molecule anchor uses B
    f = trial_forecast(0.3, 0.09, 0.1, 0.70, 0.3, 1.5); assert f["P_cf_member"] == 0.70
    assert f["P_c"] == round(float(np.median([f["P1"], 0.5, 0.70])), 4)
    # 6. write-once refuses changed content in a scratch dir
    import tempfile
    from lumen import provenance as pv
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d) / "x.json"; pv.write_once(p, {"a": 1}); pv.RERUN_DIR = pathlib.Path(d) / "r"
        pv.write_once(p, {"a": 1, "run_utc": "x"})
        try:
            pv.write_once(p, {"a": 2}); raise AssertionError("write_once did not refuse")
        except pv.ProvenanceError:
            pass
    print("SELFTEST PASS: w_hat %.3f (true 0.4), T_hat %.2f, P_b(0.3,.04,.1) %.4f < P1 %.4f" % (w_hat, T_hat, pb, P1_of(0.3, 0.04, 0.1)))


# ---------------------------------------------------------------- real freeze
def freeze(date):
    from lumen.layers import drug_effect as de
    from lumen import provenance as pv
    w, T, n_units, n_trials = fit_on_cal1()
    ln, sle = json.loads(W1_LN.read_text()), json.loads(W1_SLE.read_text())
    todo = [(k, v, "ln", "W1-LN (CRR, PREREG-W1-2026-11-21)") for k, v in ln.items()]
    todo += [(k, v, TABLE_KEY[v["table"]], "W1-SLE (%s wk %s, PREREG-W1-SLE)" % (v["endpoint"], v["primary_week"]))
             for k, v in sle.items() if k != "_meta"]
    p0s = {t: table_p0(de.table(t)) for t in {x[2] for x in todo}}
    out = {"_meta": {
        "prereg": PREREG, "item": "REC1", "round": 166,
        "event": "P(observed late log OR > 0), per active arm vs its own control; unit = trial",
        "fit": {"w_mixture": w, "T_temperature_descriptive": T, "fit_data": "experiments/lumen-v3/cal1_units.csv",
                "fit_md5": CAL1_UNITS_MD5, "n_arms": n_units, "n_trials": n_trials,
                "note": "fit only; no competitor scored on historical trials (rule 8)"},
        "van_zwet": {"weights": VZ_W.tolist(), "snr_sd": VZ_TAU.tolist(), "source": "van Zwet, Schwab & Senn 2021 Stat Med Table 1, doi 10.1002/sim.9173"},
        "p0_nominal_control_rate": {t: {"median": round(v[0], 4), "n_trials": v[1]} for t, v in p0s.items()},
        "lumen_forecast_md5": pv.FROZEN["experiments/lumen-bench/lumen_forecast.py"],
        "tables_md5": {t: pv.FROZEN[de.TABLES[t]] for t in p0s},
        "uncertainty_sources": {"P1": "between-trial tau + anticipated sampling", "P_a": "between-trial tau + anticipated sampling",
                                "P_c": "between-trial tau + anticipated sampling", "P_b": "between-trial tau + anticipated sampling + Cochrane SNR prior",
                                "B": "none (base rate)"},
        "excluded": [], "exclusion_screen": "experiments/lumen-bench/prospective/rec1/topline_screen_2026-11-28.md",
        "sources_manifest": "experiments/lumen-bench/prospective/rec1/sources/MANIFEST.tsv"}}
    for nct, rec, tab, src in todo:
        ct = json.loads((SRC / f"{nct}.json").read_text())
        assert ct.get("hasResults") is False, nct
        enr = ct["protocolSection"]["designModule"]["enrollmentInfo"]
        n_arms, ctrl = arms_from_ctgov(ct)
        pr = de.prior(rec["molecule"], rec["target"], tab, anchor=False)
        M1, sd = float(pr["prior_logOR"]), float(pr["prior_sd"]); v = sd ** 2
        if tab != "ln":
            fz = rec["M1_pooled"]
            assert abs(M1 - fz["prior_logOR"]) <= 1e-3 and abs(sd - fz["prior_sd"]) <= 1e-3, (nct, M1, sd, fz)
        pa = de.prior(rec["molecule"], rec["target"], tab, anchor=True)
        cf = (float(pa["prior_logOR"]), float(pa["prior_sd"]) ** 2) if pa["anchor"] == "molecule" else None
        B = de.forecast(rec["molecule"], rec["target"], tab)["base_rate_positive"]["value"]
        p0 = p0s[tab][0]
        n_a = enr["count"] / n_arms
        s2 = se_ant2(n_a, p0, M1)
        f = trial_forecast(M1, v, s2, B, w, T, cf)
        out[nct] = {"trial": rec["trial"], "molecule": rec["molecule"], "target": rec["target"], "source_freeze": src,
                    "table": tab, "enrollment": enr, "n_arms": n_arms, "control_arms": ctrl,
                    "scorable_note": None if ctrl else "no control arm in armGroups: unscorable unless the results name a comparator",
                    "n_per_arm_anticipated": round(n_a, 2), "se_ant": round(math.sqrt(s2), 4),
                    "M1": M1, "v": round(v, 6), "anchor_M3": pa["anchor"],
                    "M3_cf": None if cf is None else {"logOR": cf[0], "v": round(cf[1], 6)}, **f}
    out = pv.stamp(out, prereg=PREREG)
    out["_meta"]["frozen_utc"] = out.pop("run_utc"); out["_meta"]["git_head"] = out.pop("git_head")
    out["_meta"]["prereg_commit"] = out.pop("prereg_commit")
    path = OUT_DIR / f"frozen_REC1_{date}.json"
    res = pv.write_once(path, out)
    print("wrote", res, "w", w, "T", T)
    for k, r in out.items():
        if k != "_meta":
            print(k, r["table"], "M1 %.3f se_ant %.3f  P1 %.3f a %.3f b %.3f c %.3f B %.2f D %.3f" %
                  (r["M1"], r["se_ant"], r["P1"], r["P_a"], r["P_b"], r["P_c"], r["B"], r["D_legacy_latent"]))
    return path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--date", default="2026-11-28")
    a = ap.parse_args()
    selftest() if a.selftest else freeze(a.date)
