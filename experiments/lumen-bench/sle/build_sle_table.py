"""P3 data step: SRI-4 (and BICLA) responders per arm per visit for placebo-controlled non-renal SLE RCTs.
Usage: python3 -I build_sle_table.py <lupus_with_results.json.gz> <out_sri4.csv> <out_bicla.csv>
Rules: experiments/lumen-bench/sle/PREREG-P3.md. LABELS fixed in the PREREG commit."""
import gzip, json, re, sys, math, csv
# NCT -> (trial label, molecule, target, onset class)
LABELS = {
 'NCT04433585': ('LY3471851-ph2', 'LY3471851', 'IL-2R', 'other'),
 'NCT01283139': ('sifalimumab-ph2', 'sifalimumab', 'IFN-alpha', 'fast'),
 'NCT03616912': ('SLE-BRAVE-I', 'baricitinib', 'JAK', 'fast'),
 'NCT03616964': ('SLE-BRAVE-II', 'baricitinib', 'JAK', 'fast'),
 'NCT02708095': ('baricitinib-ph2', 'baricitinib', 'JAK', 'fast'),
 'NCT04925934': ('daxdilimab-ph2', 'daxdilimab', 'ILT7', 'fast'),
 'NCT04186871': ('branebrutinib-ph2', 'branebrutinib', 'BTK', 'slow'),
 'NCT01649765': ('PLUTO', 'belimumab', 'BAFF', 'slow'),
 'NCT03845517': ('brepocitinib-ph2', 'brepocitinib', 'TYK2', 'fast'),
 'NCT00424476': ('BLISS-52', 'belimumab', 'BAFF', 'slow'),
 'NCT01972568': ('ADDRESS-II', 'atacicept', 'BAFF/APRIL', 'slow'),
 'NCT02437890': ('vobarilizumab-ph2', 'vobarilizumab', 'IL-6', 'other'),
 'NCT05162586': ('WILLOW', 'enpatoran', 'TLR7/8', 'other'),
 'NCT04835441': ('acazicolcept-ph2', 'acazicolcept', 'ICOSL/CD86', 'slow'),
 'NCT02975336': ('evobrutinib-ph2', 'evobrutinib', 'BTK', 'slow'),
 'NCT02270957': ('ACCESS-abatacept-SLE', 'abatacept', 'CD80/86', 'slow'),
 'NCT02908100': ('fenebrutinib-ph2', 'fenebrutinib', 'BTK', 'slow'),
 'NCT02660944': ('RSLV-132-ph2', 'RSLV-132', 'RNA', 'other'),
 'NCT02504645': ('rigerimod-ph3', 'rigerimod', 'P140 peptide', 'other'),
 'NCT01405196': ('PF-04236921-ph2', 'PF-04236921', 'IL-6', 'other'),
 'NCT01345253': ('belimumab-NEA', 'belimumab', 'BAFF', 'slow'),
 'NCT04963296': ('ALLEGORY', 'obinutuzumab', 'CD20', 'slow'),
 'NCT01438489': ('MUSE', 'anifrolumab', 'IFNAR', 'fast'),
 'NCT02665364': ('IFN-K-ph2b', 'IFN-kinoid', 'IFN-alpha', 'fast'),
 'NCT02446899': ('TULIP-2', 'anifrolumab', 'IFNAR', 'fast'),
 'NCT02446912': ('TULIP-1', 'anifrolumab', 'IFNAR', 'fast'),
 'NCT04931563': ('anifrolumab-Asia', 'anifrolumab', 'IFNAR', 'fast'),
 'NCT02847598': ('LILAC', 'litifilimab', 'BDCA2', 'fast'),
 'NCT02265744': ('lulizumab-ph2', 'lulizumab', 'CD28', 'slow'),
 'NCT05087628': ('PRV-3279-ph2', 'PRV-3279', 'CD32B/CD79B', 'slow'),
 'NCT04294667': ('PHOENYCS-GO', 'dapirolizumab', 'CD40L', 'slow'),
 'NCT02804763': ('dapirolizumab-ph2', 'dapirolizumab', 'CD40L', 'slow'),
 'NCT01196091': ('ILLUMINATE-1', 'tabalumab', 'BAFF', 'slow'),
 'NCT01205438': ('ILLUMINATE-2', 'tabalumab', 'BAFF', 'slow'),
 'NCT03252587': ('PAISLEY', 'deucravacitinib', 'TYK2', 'fast'),
 'NCT05123586': ('LY3361237-ph2', 'LY3361237', 'BTLA', 'other'),
 'NCT01632241': ('EMBRACE', 'belimumab', 'BAFF', 'slow'),
 'NCT03371251': ('BOS161721-ph2', 'BOS161721', 'IL-21', 'other'),
 'NCT03093402': ('lenabasum-ph2', 'lenabasum', 'CB2', 'other'),
 'NCT01484496': ('BLISS-SC', 'belimumab', 'BAFF', 'slow'),
 'NCT03161483': ('iberdomide-ph2', 'iberdomide', 'cereblon', 'other'),
 'NCT01135459': ('rigerimod-ph2', 'rigerimod', 'P140 peptide', 'other'),
 'NCT03742037': ('CARE-cenerimod', 'cenerimod', 'S1P1', 'other'),
 'NCT04058028': ('rozibafusp-ph2', 'rozibafusp', 'BAFF/ICOSL', 'slow'),
 'NCT02349061': ('ustekinumab-ph2', 'ustekinumab', 'IL-12/23', 'other'),
 'NCT02533570': ('brentuximab-ph2', 'brentuximab', 'CD30', 'other'),
 'NCT03616964_': None,
 'NCT03517722': ('LOTUS', 'ustekinumab', 'IL-12/23', 'other'),
 'NCT06161116': ('MK-6194-ph2', 'MK-6194', 'IL-2R', 'other'),
 'NCT00410384': ('BLISS-76', 'belimumab', 'BAFF', 'slow'),
 'NCT04680637': ('efavaleukin-ph2', 'efavaleukin', 'IL-2R', 'other'),
}
LABELS = {k: v for k, v in LABELS.items() if v}
EXCLUDED = {'NCT03656562', 'NCT03978520', 'NCT04451772', 'NCT05306873', 'NCT03843125'}
SRI = re.compile(r'SRI ?\(?4\)?|SRI-4|SLE Responder Index|Responder Index', re.I)
BIC = re.compile(r'BICLA|BILAG.{0,20}based (Composite|Combined)', re.I)
EXC = re.compile(r'[5-8]\)|SRI-?[5-8]|steroid|OCS|corticosteroid|prednis|sustain|time to|duration|interferon|IFN|test.?high|'
                 r'test.?low|modified|mSRI|subgroup|stratum|composite response|with no|without|and OCG|dose ≤|dose <=', re.I)
DROPG = re.compile(r'open-label|extension|\bLTE\b|long-term|safety follow|Part B|Part C|\(MEE\)|Cohort A|Phase 1b|pSS|'
                   r'Sj[oö]gren|Rheumatoid|\(RA\)|Main Study\) /', re.I)
DOSE = re.compile(r'\d+(\.\d+)? ?(mg|mcg|µg|mg/kg)', re.I)
def weeks(txt):
    m = re.search(r'(week|day|month)s?\s*(\d+)', txt or '', re.I)
    if not m: return None
    u, v = m.group(1).lower(), int(m.group(2))
    return v if u == 'week' else round(v / 7) if u == 'day' else round(v * 52 / 12)
def half_up(x): return int(math.floor(x + 0.5))
D = {s['protocolSection']['identificationModule']['nctId']: s for s in json.load(gzip.open(sys.argv[1]))}
def build(endpoint, extra_exc=None):
    cand = []
    for nct, (lab, mol, tgt, cls) in LABELS.items():
        if nct in EXCLUDED or nct not in D: continue
        for i, m in enumerate(D[nct].get('resultsSection', {}).get('outcomeMeasuresModule', {}).get('outcomeMeasures', [])):
            ttl = m['title']
            if not endpoint.search(ttl): continue
            if EXC.search(ttl) and not (endpoint is BIC and re.search(r'BICLA', ttl) and not EXC.search(ttl.replace('BICLA', ''))): continue
            if endpoint is SRI and BIC.search(ttl): continue
            if m.get('paramType') not in ('COUNT_OF_PARTICIPANTS', 'NUMBER'): continue
            pct = 'percent' in (m.get('unitOfMeasure') or '').lower() or m.get('paramType') == 'NUMBER' and 'participants' not in (m.get('unitOfMeasure') or '').lower()
            g = {x['id']: x['title'] for x in m.get('groups', [])}
            den = {c['groupId']: int(c['value']) for d in m.get('denoms', []) for c in d['counts']}
            for c in m.get('classes', []):
                ct = c.get('title', '')
                if ct and re.search(r'non-?responder|did not|^no\b|not achiev', ct, re.I): continue
                wk = weeks(ct) or weeks(ttl) or weeks(m.get('timeFrame', ''))
                if wk is None: continue
                for cat in c['categories']:
                    if cat.get('title') and re.search(r'non-?responder|^no$|did not|not achiev', cat['title'], re.I): continue
                    for me in cat['measurements']:
                        v = me.get('value'); gid = me['groupId']
                        if v in (None, 'NA') or gid not in den or den[gid] == 0: continue
                        N = den[gid]; r = half_up(float(v) * N / 100) if pct else int(float(v))
                        cand.append((nct, g[gid], wk, r, N, i, m['type'], ttl))
    best = {}
    for nct, gt, wk, r, N, i, ty, ttl in cand:
        if DROPG.search(gt): continue
        is_ctrl = bool(re.search(r'placebo', gt, re.I)) and not DOSE.search(gt)
        if is_ctrl and re.search(r'followed by', gt, re.I) and wk > 24: continue
        arm = 'ctrl' if is_ctrl else re.sub(r'\s+', ' ', gt)[:60]
        k = (nct, arm, wk); sc = (N, ty == 'PRIMARY', -i)
        if k not in best or sc > best[k][0]: best[k] = (sc, r, N, i, ttl)
    maxN = {}
    for (nct, a, w), v in best.items(): maxN[(nct, a)] = max(maxN.get((nct, a), 0), v[2])
    rows = []
    for (nct, a, w), v in sorted(best.items()):
        if a == 'ctrl': continue
        c = best.get((nct, 'ctrl', w))
        if c is None or v[2] < 0.5 * maxN[(nct, a)] or c[2] < 0.5 * maxN[(nct, 'ctrl')]: continue
        lab, mol, tgt, cls = LABELS[nct]
        rows.append(dict(trial=lab, nct=nct, arm=a, molecule=mol, target=tgt, cls=cls, week=w, r_drug=v[1], n_drug=v[2],
                         r_ctrl=c[1], n_ctrl=c[2], measure_drug=v[3], measure_ctrl=c[3], title=v[4][:100]))
    return rows
for ep, out in [(SRI, sys.argv[2]), (BIC, sys.argv[3])]:
    rows = build(ep)
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(out, len(rows), 'rows,', len({(r['trial'], r['arm']) for r in rows}), 'arms,', len({r['trial'] for r in rows}), 'trials')
