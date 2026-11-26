"""RT-20261126-1 run 1 (round 146): catalogue of unreported non-renal SLE RCTs for a prospective freeze.
Registry METADATA only: no outcome numbers are read (studies with hasResults=True are dropped before any field is used).
Input: raw/ctgov_sle_p1.json (fetched with scripts/provenance/fetch.sh, see MANIFEST.tsv).
Output: ctgov_screen.csv (every study + the first rule it failed) and candidates_w1sle.csv (passed registry rules 1-7).
Rule 8 (no public topline before the PREREG / before June 2026) is a manual screen recorded in topline_screen.csv.
"""
import json, re, csv, os
D = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(D, '..', '..')

# trials whose outcomes LUMEN has already tabulated or looked at (not prospective)
seen = set()
for f in ['resp_table.csv', 'sle/sri4_table.csv', 'sle/bicla_table.csv', 'regdocs/rt4_rows.csv',
          'newtrials/candidates_frozen.csv', 'newtrials/ctgov_candidates.csv']:
    p = os.path.join(B, f)
    if os.path.exists(p):
        seen |= set(re.findall(r'NCT\d{8}', open(p).read()))
w1_ln = set(json.load(open(os.path.join(B, 'prospective/w1/w1_frozen_2026-11-21.json'))).keys())

RESP = re.compile(r'respon|remission|\bSRI\b|SRI-?\d|BICLA|LLDAS|DORIS|CLASI-?50|\bCRR\b|\bPRR\b', re.I)
CTRL = re.compile(r'placebo|standard (of )?care|\bSoC\b|background', re.I)

def weeks(tf):
    tf = tf.lower(); out = []
    for num, unit in re.findall(r'(\d+(?:\.\d+)?)\s*-?\s*(weeks?|wks?|months?|mos?|years?|yrs?|days?)', tf):
        x = float(num)
        out.append(x if unit.startswith('w') else x * 52 / 12 if unit.startswith('mo') else x * 52 if unit.startswith('y') else x / 7)
    out += [float(n) for n in re.findall(r'\bweek\s*(\d+)', tf)]
    out += [float(n) * 52 / 12 for n in re.findall(r'\bmonth\s*(\d+)', tf)]
    return sorted(set(round(w) for w in out))

screen, cands = [], []
for s in json.load(open(os.path.join(D, 'raw', 'ctgov_sle_p1.json')))['studies']:
    p = s['protocolSection']; idm = p['identificationModule']; nct = idm['nctId']
    title = idm.get('briefTitle', ''); conds = ';'.join(p.get('conditionsModule', {}).get('conditions', []))
    dm = p.get('designModule', {}); st = p.get('statusModule', {})
    arms = p.get('armsInterventionsModule', {}).get('armGroups', [])
    def fail(rule):
        screen.append(dict(nct=nct, title=title, failed_rule=rule)); return True
    if s.get('hasResults'): fail('1 results posted on CT.gov'); continue
    if nct in seen: fail('2 already tabulated/seen by LUMEN'); continue
    if nct in w1_ln: fail('2b already in W1 LN freeze'); continue
    if re.search(r'nephrit|kidney|renal', conds + ' ' + title, re.I): fail('3 renal (covered by W1 LN)'); continue
    if re.search(r'cutaneous|discoid|\bCLE\b|skin', title, re.I) and not re.search(r'systemic', title, re.I):
        fail('3b cutaneous-only'); continue
    if dm.get('designInfo', {}).get('allocation') != 'RANDOMIZED': fail('4 not randomised'); continue
    ctrl = [a for a in arms if a.get('type') in ('PLACEBO_COMPARATOR', 'NO_INTERVENTION', 'SHAM_COMPARATOR')
            or CTRL.search(a.get('label', '') + ' ' + a.get('description', '')[:300]) and a.get('type') != 'EXPERIMENTAL']
    if not ctrl: fail('5 no placebo/SoC arm'); continue
    n = (dm.get('enrollmentInfo', {}) or {}).get('count', 0) or 0
    if n < 20: fail('6 enrolled < 20'); continue
    prim = p.get('outcomesModule', {}).get('primaryOutcomes', [])
    resp = [o for o in prim if RESP.search(o.get('measure', '') + ' ' + o.get('description', '')[:300])]
    if not resp: fail('7 primary endpoint not a responder proportion (not scorable as log OR)'); continue
    o = resp[0]
    refs = p.get('referencesModule', {}).get('references', [])
    r = dict(nct=nct, title=title, phase='/'.join(dm.get('phases', [])), status=st.get('overallStatus', ''),
             primary_completion=st.get('primaryCompletionDateStruct', {}).get('date', ''),
             completion=st.get('completionDateStruct', {}).get('date', ''), enrolled=n,
             sponsor=p.get('sponsorCollaboratorsModule', {}).get('leadSponsor', {}).get('name', ''),
             interventions=';'.join(sorted({i.get('name', '') for i in p.get('armsInterventionsModule', {}).get('interventions', [])})),
             n_arms=len(arms), ctrl_arms=';'.join(a.get('label', '') for a in ctrl),
             primary_measure=o.get('measure', '')[:200], primary_timeframe=o.get('timeFrame', '')[:120],
             primary_weeks=';'.join(str(w) for w in weeks(o.get('timeFrame', ''))),
             other_responder_weeks=';'.join(str(w) for w in sorted({w for x in p.get('outcomesModule', {}).get('secondaryOutcomes', [])
                                             if RESP.search(x.get('measure', '')) for w in weeks(x.get('timeFrame', ''))})),
             registry_result_refs=';'.join(x.get('pmid', '') for x in refs if x.get('type') == 'RESULT'))
    cands.append(r); screen.append(dict(nct=nct, title=title, failed_rule='passed 1-7'))

with open(os.path.join(D, 'ctgov_screen.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['nct', 'title', 'failed_rule']); w.writeheader(); w.writerows(sorted(screen, key=lambda r: r['nct']))
with open(os.path.join(D, 'candidates_w1sle.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(cands[0].keys())); w.writeheader(); w.writerows(sorted(cands, key=lambda r: r['primary_completion']))
from collections import Counter
print(Counter(r['failed_rule'] for r in screen)); print('candidates:', len(cands))
