"""lumen CLI.
  lumen forecast MOLECULE TARGET [--table ln|sri4|bicla] [--early WEEK R_DRUG N_DRUG R_CTRL N_CTRL] [--anchor] [--holdout TRIAL]
  lumen status        which modules ship, at which evidence label and n
  lumen verify        md5 of every frozen input
Research and simulation only. Not clinical guidance.
"""
from __future__ import annotations
import argparse, json, sys
from .layers import drug_effect
from . import registry, provenance


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lumen")
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("forecast")
    f.add_argument("molecule"); f.add_argument("target")
    f.add_argument("--table", default="ln", choices=sorted(drug_effect.TABLES))
    f.add_argument("--holdout"); f.add_argument("--anchor", action="store_true")
    f.add_argument("--early", nargs=5, type=float, metavar=("WEEK", "R_DRUG", "N_DRUG", "R_CTRL", "N_CTRL"))
    f.add_argument("--use-early-update", action="store_true")
    sub.add_parser("status"); sub.add_parser("verify")
    x = ap.parse_args(argv)
    if x.cmd == "forecast":
        out = drug_effect.forecast(x.molecule, x.target, x.table, x.early, x.holdout, x.anchor, x.use_early_update)
    elif x.cmd == "status":
        out = registry.as_table()
    else:
        out = provenance.verify_all()
    json.dump(out, sys.stdout, indent=1); print()
    return 0 if x.cmd != "verify" or all(out.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
