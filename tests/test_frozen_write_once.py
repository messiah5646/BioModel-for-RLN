"""Prospective freezes are write-once (reviewer R-20261127-2, round 158)."""
import json, pathlib
import pytest
from lumen import provenance

P = provenance.BENCH / "prospective"


def test_every_freeze_listed_and_unchanged():
    want = provenance.frozen_md5_table()
    found = {str(p.relative_to(provenance.REPO)) for p in P.rglob("*frozen*.json")}
    assert found == set(want), f"unlisted or missing freezes: {found ^ set(want)} (append new ones to FROZEN-MD5.tsv)"
    for rel, h in want.items():
        assert provenance.md5(provenance.REPO / rel) == h, rel


def test_write_once_never_overwrites(tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "RERUN_DIR", tmp_path / "reruns")
    f = tmp_path / "frozen_x.json"
    provenance.write_frozen(f, {"a": 1, "run_utc": "t0"})
    before = f.read_bytes()
    out = provenance.write_frozen(f, {"a": 1, "run_utc": "t1", "git_head": "abc"})   # stamps differ only: OK
    assert f.read_bytes() == before and out.name == "frozen_x.rerun.json"
    with pytest.raises(provenance.ProvenanceError):
        provenance.write_frozen(f, {"a": 2, "run_utc": "t2"})                       # content differs: refuse
    assert f.read_bytes() == before
