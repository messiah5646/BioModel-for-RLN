"""RT-20261124-1 run 1: catalogue (registry METADATA only; no outcome numbers are read).
Selects randomised LN/SLE trials with NO results on ClinicalTrials.gov, not among the 60 trials
already tabulated, whose REGISTERED responder outcomes include an early (12-26 wk) and a late
(48-104 wk) time point, and which list a RESULT-type publication (PMID). Output: candidates CSV."""
import json, re, csv, sys, os
D = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(D, '..')
known = set()
for f in ['resp_table.csv', 'sle/sri4_table.csv', 'sle/bicla_table.csv']:
    for r in csv.DictReader(open(os.path.join(B, f))): known.add(r['nct'])
RESP = re.compile(r'respon|remission|\bSRI\b|SRI-?4|BICLA|\bCRR\b|\bPRR\b|\bCR\b|\bPR\b', re.I)
def weeks(tf):
    tf = tf.lower(); out = []
    for num, unit in re.findall(r'(\d+(?:\.\d+)?)\s*-?\s*(weeks?|wks?|months?|mos?|years?|yrs?|days?)', tf):
        x = float(num)
        if unit.startswith('w'): out.append(x)
        elif unit.startswith('mo'): out.append(x * 52 / 12)
        elif unit.startswith('y'): out.append(x * 52)
        elif unit.startswith('d'): out.append(x / 7)
    for num in re.findall(r'\bweek\s*(\d+)', tf): out.append(float(num))
    for num in re.findall(r'\bmonth\s*(\d+)', tf): out.append(float(num) * 52 / 12)
    return out
rows = {}
for fn in sorted(os.listdir(os.path.join(D, 'raw'))):
    if not fn.endswith('.json'): continue
    for s in json.load(open(os.path.join(D, 'raw', fn)))['studies']:
        p = s['protocolSection']; nct = p['identificationModule']['nctId']
        if nct in known or nct in rows: continue
        ph = p.get('designModule', {}).get('phases', [])
        n = p.get('designModule', {}).get('enrollmentInfo', {}).get('count', 0) or 0
        om = p.get('outcomesModule', {})
        tps = []
        for o in om.get('primaryOutcomes', []) + om.get('secondaryOutcomes', []):
            if RESP.search(o.get('measure', '') + ' ' + o.get('description', '')[:200]):
                tps += weeks(o.get('timeFrame', ''))
        early = any(12 <= w <= 26.5 for w in tps); late = any(47.5 <= w <= 104.5 for w in tps)
        refs = p.get('referencesModule', {}).get('references', [])
        pm = [r.get('pmid') for r in refs if r.get('type') == 'RESULT' and r.get('pmid')]
        conds = ';'.join(p.get('conditionsModule', {}).get('conditions', []))
        ln = bool(re.search(r'nephrit', conds + p['identificationModule']['briefTitle'], re.I))
        rows[nct] = dict(nct=nct, indication='LN' if ln else 'SLE', title=p['identificationModule']['briefTitle'],
            phase='/'.join(ph), enrolled=n, sponsor=p.get('sponsorCollaboratorsModule', {}).get('leadSponsor', {}).get('name', ''),
            countries=';'.join(sorted({l.get('country', '') for l in p.get('contactsLocationsModule', {}).get('locations', [])})),
            responder_tp_weeks=';'.join(str(round(w)) for w in sorted(set(tps))), early=early, late=late,
            result_pmids=';'.join(pm))
keep = [r for r in rows.values() if r['early'] and r['late'] and r['enrolled'] >= 20
        and any(x in r['phase'] for x in ('PHASE2', 'PHASE3', 'PHASE4'))]
print('not-in-60, no CT.gov results:', len(rows), '| early+late responder timepoints, phase>=2, n>=20:', len(keep),
      '| with RESULT PMID:', sum(1 for r in keep if r['result_pmids']))
with open(os.path.join(D, 'ctgov_candidates.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(keep[0].keys())); w.writeheader()
    for r in sorted(keep, key=lambda r: (r['indication'], r['nct'])): w.writerow(r)
