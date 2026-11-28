"""W3 run 1 (round 169): catalogue of unreported SLE + LN RCTs beyond the 21 W1/REC1 trials, for the rolling freeze.
Registry METADATA only: studies with hasResults=True are dropped before any other field is used; no outcome is read.
Input: raw/ctgov_lupus_open.json (RECRUITING|NOT_YET_RECRUITING|ENROLLING_BY_INVITATION) and raw/ctgov_lupus_closed.json
(ACTIVE_NOT_RECRUITING|COMPLETED), fetched via scripts/provenance/fetch.sh (MANIFEST.tsv). Rules mirror
w1-sle/catalogue_w1sle.py, extended to renal trials and open statuses.
Output: w3_screen.csv (first failed rule per study), w3_candidates.csv (passed registry rules), accrual table on stdout.
"""
import json, re, csv, os
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(D, '..', '..')
seen = set()
for f in ['resp_table.csv', 'sle/sri4_table.csv', 'sle/bicla_table.csv', 'regdocs/rt4_rows.csv',
          'newtrials/candidates_frozen.csv', 'newtrials/ctgov_candidates.csv']:
    p = os.path.join(B, f)
    if os.path.exists(p):
        seen |= set(re.findall(r'NCT\d{8}', open(p).read()))
w21 = set(json.load(open(os.path.join(B, 'prospective/w1/w1_frozen_2026-11-21.json'))).keys())
w21 |= set(re.findall(r'NCT\d{8}', open(os.path.join(B, 'prospective/rec1/frozen_REC1_2026-11-28.json')).read()))
# W1-SLE run-1 manual exclusions (public topline before 2026-11-26) stay excluded
w1sle_top = set(re.findall(r'NCT\d{8}', open(os.path.join(B, 'prospective/w1-sle/topline_screen.csv')).read()))

RESP = re.compile(r'respon|remission|\bSRI\b|SRI-?\d|BICLA|\bCRR\b|\bPRR\b|\bORR\b', re.I)
NONDA = re.compile(r'vaccin|LLDAS|DORIS|CLASI|pharmacokinet|safety|adverse|pregnan|fatigue|cardiovasc|bone|thromb', re.I)
CTRL = re.compile(r'placebo|standard (of )?care|\bSoC\b|background|conventional', re.I)
SLE = re.compile(r'systemic lupus|\bSLE\b|lupus nephritis|\bLN\b|nephritis|lupus', re.I)

def weeks(tf):
    tf = tf.lower(); out = []
    for num, unit in re.findall(r'(\d+(?:\.\d+)?)\s*-?\s*(weeks?|wks?|months?|mos?|years?|yrs?|days?)', tf):
        x = float(num)
        out.append(x if unit.startswith('w') else x * 52 / 12 if unit.startswith('mo') else x * 52 if unit.startswith('y') else x / 7)
    out += [float(n) for n in re.findall(r'\bweek\s*(\d+)', tf)]
    out += [float(n) * 52 / 12 for n in re.findall(r'\bmonth\s*(\d+)', tf)]
    return sorted(set(round(w) for w in out))

screen, cands = [], []
studies = {}
for t in ['open', 'closed']:
    for s in json.load(open(os.path.join(D, 'raw', f'ctgov_lupus_{t}.json')))['studies']:
        studies[s['protocolSection']['identificationModule']['nctId']] = s
for nct, s in sorted(studies.items()):
    p = s['protocolSection']; idm = p['identificationModule']
    title = idm.get('briefTitle', ''); conds = ';'.join(p.get('conditionsModule', {}).get('conditions', []))
    dm = p.get('designModule', {}); st = p.get('statusModule', {})
    arms = p.get('armsInterventionsModule', {}).get('armGroups', [])
    def fail(rule): screen.append(dict(nct=nct, title=title, failed_rule=rule))
    if s.get('hasResults'): fail('1 results posted on CT.gov'); continue
    if nct in seen: fail('2 already tabulated/seen by LUMEN'); continue
    if nct in w21: fail('2b already in W1/REC1 21-trial freeze'); continue
    if nct in w1sle_top: fail('2c W1-SLE manual exclusion (topline/endpoint) 2026-11-26'); continue
    if not re.search(r'systemic lupus|\bSLE\b|lupus nephritis', conds + ' ' + title, re.I): fail('3 not SLE/LN'); continue
    if re.search(r'cutaneous|discoid|\bCLE\b|skin', title, re.I) and not re.search(r'systemic', title, re.I):
        fail('3b cutaneous-only'); continue
    if dm.get('studyType') != 'INTERVENTIONAL' or dm.get('designInfo', {}).get('allocation') != 'RANDOMIZED':
        fail('4 not randomised interventional'); continue
    ctrl = [a for a in arms if a.get('type') in ('PLACEBO_COMPARATOR', 'NO_INTERVENTION', 'SHAM_COMPARATOR')
            or CTRL.search(a.get('label', '') + ' ' + a.get('description', '')[:300]) and a.get('type') != 'EXPERIMENTAL']
    if not ctrl: fail('5 no placebo/SoC arm'); continue
    n = (dm.get('enrollmentInfo', {}) or {}).get('count', 0) or 0
    if n < 20: fail('6 enrolled < 20'); continue
    prim = p.get('outcomesModule', {}).get('primaryOutcomes', [])
    resp = [o for o in prim if RESP.search(o.get('measure', '') + ' ' + o.get('description', '')[:300])]
    if not resp: fail('7 primary endpoint not a responder proportion'); continue
    o = resp[0]
    renal = bool(re.search(r'nephrit|kidney|renal', conds + ' ' + title, re.I))
    cands.append(dict(nct=nct, title=title, renal=int(renal), phase='/'.join(dm.get('phases', [])),
        status=st.get('overallStatus', ''), first_posted=st.get('studyFirstPostDateStruct', {}).get('date', ''),
        primary_completion=st.get('primaryCompletionDateStruct', {}).get('date', ''), enrolled=n,
        sponsor=p.get('sponsorCollaboratorsModule', {}).get('leadSponsor', {}).get('name', ''),
        interventions=';'.join(sorted({i.get('name', '') for i in p.get('armsInterventionsModule', {}).get('interventions', [])})),
        n_arms=len(arms), ctrl_arms=';'.join(a.get('label', '') for a in ctrl),
        primary_measure=o.get('measure', '')[:200], primary_timeframe=o.get('timeFrame', '')[:120],
        primary_weeks=';'.join(str(w) for w in weeks(o.get('timeFrame', ''))),
        nonDA_flag=int(bool(NONDA.search(o.get('measure', ''))))))
    screen.append(dict(nct=nct, title=title, failed_rule='passed 1-7'))

with open(os.path.join(D, 'w3_screen.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['nct', 'title', 'failed_rule']); w.writeheader(); w.writerows(screen)
with open(os.path.join(D, 'w3_candidates.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(cands[0].keys())); w.writeheader(); w.writerows(sorted(cands, key=lambda r: r['primary_completion']))
print(Counter(r['failed_rule'] for r in screen)); print('candidates:', len(cands), 'renal:', sum(c['renal'] for c in cands))
print('first-posted year:', sorted(Counter(c['first_posted'][:4] for c in cands).items()))
print('status:', Counter(c['status'] for c in cands))
