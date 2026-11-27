"""L0 provenance: repo paths, md5 checks for frozen inputs, and run stamps (rule 12b)."""
from __future__ import annotations
import datetime, hashlib, json, pathlib, subprocess

REPO = pathlib.Path(__file__).resolve().parents[1]
BENCH = REPO / "experiments" / "lumen-bench"

# Inputs whose bytes are part of published / frozen results. Changing one is a model change (rule 12b amendment).
FROZEN = {
    "experiments/lumen-bench/lumen_forecast.py": "2c11b9a5aba5c047c0817140c7d5cbe6",
    "experiments/lumen-bench/resp_table.csv": "3ce67bdb51c6597bce6666d914618778",
    "experiments/lumen-bench/resp_table_plus_pub.csv": "9631ac1b0638f5cfd17ffa017a457b8a",
    "experiments/lumen-bench/sle/sri4_table.csv": "b37dd7cd72c9e912ba079088532c3765",
    "experiments/lumen-bench/sle/bicla_table.csv": "99eb279e207302973fef6899f0e86722",
}

class ProvenanceError(RuntimeError):
    pass

def md5(path) -> str:
    return hashlib.md5(pathlib.Path(path).read_bytes()).hexdigest()

def verify(rel: str) -> pathlib.Path:
    """Return the absolute path of a frozen input after checking its md5; raise if it has drifted."""
    p = REPO / rel
    got = md5(p)
    if got != FROZEN[rel]:
        raise ProvenanceError(f"{rel}: md5 {got} != frozen {FROZEN[rel]} (a model change needs a PREREG amendment)")
    return p

def verify_all() -> dict:
    return {rel: md5(REPO / rel) == want for rel, want in FROZEN.items()}

def _git(*a) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(REPO), *a], text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "NA"

def stamp(result: dict, prereg: str | None = None) -> dict:
    """Same fields as scripts/provenance/stamp.py: run_utc, git_head, prereg_commit."""
    out = dict(result)
    out["run_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    out["git_head"] = _git("rev-parse", "--short", "HEAD")
    if prereg:
        out["prereg_commit"] = _git("log", "-1", "--format=%h %cI", "--", prereg)
    return out


# --- Write-once prospective freezes (reviewer R-20261127-2, round 158) -----------------------------------------
# A frozen prospective forecast is a public prediction: once committed its bytes must never change. Re-running a
# script that made one writes <name>.rerun.json to a scratch dir and checks the forecast content is identical.
STAMP_KEYS = frozenset({"run_utc", "git_head", "prereg_commit", "frozen_utc"})
FROZEN_MD5_TSV = BENCH / "prospective" / "FROZEN-MD5.tsv"
RERUN_DIR = REPO / ".reruns"   # git-ignored scratch


def _strip_stamps(x):
    if isinstance(x, dict):
        return {k: _strip_stamps(v) for k, v in x.items() if k not in STAMP_KEYS}
    if isinstance(x, list):
        return [_strip_stamps(v) for v in x]
    return x


def write_once(path, obj, indent=1, default=None) -> pathlib.Path:
    """Write a prospective freeze (or a committed result JSON) once. If `path` exists: never overwrite; write
    <stem>.rerun.json under RERUN_DIR and raise ProvenanceError unless the content (stamp keys removed) equals the
    committed file."""
    path = pathlib.Path(path)
    text = json.dumps(obj, indent=indent, default=default)
    if not path.exists():
        path.write_text(text)
        return path
    RERUN_DIR.mkdir(parents=True, exist_ok=True)
    rerun = RERUN_DIR / (path.stem + ".rerun.json")
    rerun.write_text(text)
    if _strip_stamps(json.loads(text)) != _strip_stamps(json.loads(path.read_text())):
        raise ProvenanceError(f"{path.name}: re-run forecasts differ from the committed freeze (see {rerun})")
    return rerun


write_frozen = write_once


def write_once_text(path, text: str) -> pathlib.Path:
    """Same rule for text outputs (prediction CSVs): byte-identical re-runs are a no-op, differences raise."""
    path = pathlib.Path(path)
    if not path.exists():
        path.write_text(text)
        return path
    if path.read_text() == text:
        return path
    RERUN_DIR.mkdir(parents=True, exist_ok=True)
    rerun = RERUN_DIR / (path.stem + ".rerun" + path.suffix)
    rerun.write_text(text)
    raise ProvenanceError(f"{path.name}: re-run output differs from the committed file (see {rerun})")


def frozen_md5_table() -> dict:
    rows = [l.split("\t") for l in FROZEN_MD5_TSV.read_text().splitlines() if l and not l.startswith("#")]
    return {rel: h for rel, h in rows}
