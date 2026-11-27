import numpy as np, pandas as pd, pytest
from lumen import gate, baselines as bl, registry
from lumen.contracts import Observable, ArmCount, ContractError, combinable, check_visits
from lumen.provenance import BENCH


def test_gate_fail_branch_rule_2b():
    r = gate.verdict(0.259, {"M0": 0.367, "M1": 0.311, "C0": 0.146}, "M1", 13, 7)
    assert r.verdict == "FAIL"           # RT-1 numbers: shrunk update vs carry-forward


def test_gate_pass_and_underpowered():
    assert gate.verdict(0.08, {"M0": 0.17, "M1": 0.11}, "M1", 30, 20).verdict == "PASS"
    assert gate.verdict(0.08, {"M0": 0.17, "M1": 0.11}, "M1", 5, 3).verdict == "UNDERPOWERED"
    assert gate.verdict(0.105, {"M0": 0.17, "M1": 0.11}, "M1", 30, 20).verdict == "INCONCLUSIVE"


def test_scores_basic():
    obs, se = np.array([0.3, -0.1]), np.array([0.2, 0.3])
    assert bl.iv_mse(obs, se, obs) == 0
    assert bl.coverage(obs, se, obs, np.zeros(2)) == 1.0
    assert bl.log_score(obs, se, obs, np.full(2, 0.5)) > bl.log_score(obs, se, obs, np.full(2, 0.1))


def test_contracts():
    a = Observable(1.2, "absolute", "se", 157, "ITT", "g/g", "x", spread=0.1)
    assert abs(a.sd() - 0.1 * 157 ** .5) < 1e-12
    with pytest.raises(ContractError):
        Observable(99999.0, "absolute", "sd", 10, "ITT", "g/g", "x")
    b = Observable(0.8, "ratio_to_baseline", "none", 50, "ITT", "g/g", "x")
    with pytest.raises(ContractError):
        combinable(a, b)
    with pytest.raises(ContractError):
        ArmCount("T", "a", 52, 11, 10, 1, 10)
    assert len(check_visits([ArmCount("T", "a", 24, 5, 20, 3, 20), ArmCount("T", "a", 52, 5, 20, 3, 20)])) == 1
    with pytest.raises(ContractError):
        check_visits([ArmCount("T", "a", 24, 5, 20, 3, 20), ArmCount("T", "a", 52, 5, 25, 3, 20)])


def test_armcount_matches_frozen_loader():
    R = pd.read_csv(BENCH / "resp_table.csv").iloc[0]
    lor, se = ArmCount(R.trial, R.arm, R.week, R.r_drug, R.n_drug, R.r_ctrl, R.n_ctrl).log_or()
    from lumen.layers.drug_effect import table
    T = table("ln").iloc[0]
    assert abs(lor - T.lor) < 1e-12 and abs(se - T.se) < 1e-12


def test_registry_only_pooled_ships():
    assert [m.name for m in registry.shipped()] == ["drug_effect.pooled_prior"]
