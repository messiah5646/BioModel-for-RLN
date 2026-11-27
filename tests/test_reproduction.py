"""v3 must reproduce every frozen forecast exactly: it wraps lumen_forecast.py, it does not re-implement it."""
import json, subprocess, sys
import pytest
from lumen import provenance
from lumen.layers import drug_effect as de

B = provenance.BENCH


def test_frozen_inputs_unchanged():
    assert all(provenance.verify_all().values())


def test_r90_frozen():
    F = json.load(open(B / "prospective/r90_frozen.json"))
    ins = {"NOBILITY_obinutuzumab_wk104": ("obinutuzumab", "CD20", (52, 22, 63, 14, 62), "NOBILITY"),
           "TULIP-LN1_anifro_int_wk104": ("anifrolumab", "IFNAR", (52, 23, 51, 15, 49), "TULIP-LN1"),
           "TULIP-LN1_anifro_basic_wk104": ("anifrolumab", "IFNAR", (52, 7, 45, 15, 49), "TULIP-LN1")}
    for k, (m, t, e, h) in ins.items():
        f = de.core().forecast(m, t, e, R=de.table("ln"), holdout=h)   # anchor=True, as frozen
        assert (f["post_logOR"], f["post_sd"]) == (F[k]["post_logOR"], F[k]["post_sd"])


def test_w1_sle_frozen():
    F = json.load(open(B / "prospective/w1-sle/frozen_2026-11-26.json"))
    tab = {"SRI-4": "sri4", "BICLA": "bicla"}
    for nct, v in F.items():
        if nct == "_meta":
            continue
        m1 = de.prior(v["molecule"], v["target"], tab[v["endpoint"]], anchor=False)
        m3 = de.prior(v["molecule"], v["target"], tab[v["endpoint"]], anchor=True)
        for got, want in ((m1, v["M1_pooled"]), (m3, v["M3_anchored"])):
            for k in ("anchor", "prior_logOR", "prior_sd", "OR", "OR_95", "P_benefit"):
                assert got[k] == want[k], (nct, k)


@pytest.mark.parametrize("args", [["voclosporin", "calcineurin"],
                                  ["obinutuzumab", "CD20", "--early", "52", "22", "63", "14", "62"]])
def test_cli_matches_legacy_script(args):
    legacy = json.loads(subprocess.check_output([sys.executable, "lumen_forecast.py", *args], cwd=B, text=True))
    early = None
    if "--early" in args:
        early = [float(a) for a in args[args.index("--early") + 1:]]
    v3 = de.forecast(args[0], args[1], "ln", early)
    v3.pop("layer"); v3.pop("early_note", None); legacy.pop("early_note", None)   # note text is abridged in v3
    assert v3 == legacy


def test_early_has_no_headline():
    out = de.forecast("obinutuzumab", "CD20", "ln", [52, 22, 63, 14, 62])
    assert out["headline"].startswith("NONE")
    assert set(out["early_side_by_side"]) == {"M0_null", "M1_pooled_prior", "C0_carry_forward"}
