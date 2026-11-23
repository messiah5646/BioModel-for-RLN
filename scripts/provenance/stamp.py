"""Rule 12b helper: add run provenance to a result dict before writing it to JSON.

    from stamp import stamp; json.dump(stamp(result, prereg='experiments/x/PREREG.md'), f)

Adds run_utc (time the analysis ran), git_head, and the commit time of the PREREG file,
so a reviewer can check PREREG commit < data fetch_utc < run_utc.
"""
import datetime, subprocess

def _git(*a):
    try:
        top = subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], text=True,
                                      stderr=subprocess.DEVNULL).strip()
        return subprocess.check_output(['git', '-C', top, *a], text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return 'NA'

def stamp(result, prereg=None):
    result = dict(result)
    result['run_utc'] = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    result['git_head'] = _git('rev-parse', '--short', 'HEAD')
    if prereg:  # path relative to the repo root
        result['prereg_commit'] = _git('log', '-1', '--format=%h %cI', '--', prereg)
    return result
