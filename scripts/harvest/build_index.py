"""Build per-sample index + per-series label flags from geo_headers output.
Usage: python3 -I build_index.py <hdr_dir> <catalogue.csv> <out_dir>"""
import sys, glob, re, os, pandas as pd
H, cat, O = sys.argv[1:4]
parts = []
for f in glob.glob(f'{H}/*.tsv'):
    try: g = pd.read_csv(f, sep='\t', dtype=str)
    except pd.errors.EmptyDataError: continue
    chs = [c for c in g.columns if c.startswith('ch:')]
    g['characteristics'] = ['; '.join(f'{c[3:]}={v}' for c, v in zip(chs, r) if isinstance(v, str) and v) for r in g[chs].itertuples(index=False)] if chs else ''
    g['char_keys'] = ' '.join(c[3:] for c in chs)
    pid = [c for c in chs if re.search(r'^ch:(patient|donor|subject|individual)', c)]
    g['n_pat'] = g[pid[0]].nunique() if pid else None
    parts.append(g[[c for c in ['series','geo_accession','platform_id','organism_ch1','source_name_ch1','title','characteristics','char_keys','n_pat'] if c in g]])
S = pd.concat(parts, ignore_index=True)
idx = S[['series','geo_accession','platform_id','organism_ch1','source_name_ch1','title','characteristics']]
idx.to_csv(f'{O}/geo_sample_index.tsv.gz', sep='\t', index=False)
FLAGS = {
 'nephritis': r'nephrit|\bLN\b|renal|kidney|glomerul|isn/rps|class (iii|iv|v)\b|biopsy class',
 'sledai': r'sledai|bilag|disease activity|\bsdi\b|\bpga\b',
 'response': r'respon|remission|refractor|\bCRR?\b|\bPRR?\b|outcome',
 'proteinuria': r'protein|upcr|uacr|\bpcr\b|creatinin|egfr|\bgfr\b',
 'histology': r'activity index|chronicity|\bAI\b|\bCI\b|isn|rps|crescent|fibros',
 'timepoint': r'visit|week|month|\bday\b|time ?point|baseline|follow|pre-?treat|post-?treat|days_from',
 'treatment': r'treat|drug|therapy|mmf|mycophenol|cyclophosph|rituximab|belimumab|azathioprine|steroid|prednis|hydroxychloroquine|anifrolumab|car-?t',
 'autoantibody': r'dsdna|anti-?sm|complement|\bc3\b|\bc4\b|\bana\b',
}
cat_df = pd.read_csv(cat)
rows = []
for s, g in S.groupby('series'):
    keys = str(g.char_keys.iloc[0])
    vals = ' '.join(g['characteristics'].astype(str).head(400)).lower()
    text = (keys + ' ' + vals)
    r = {'series': s, 'n_samples': len(g), 'organism': ';'.join(sorted(g.organism_ch1.dropna().unique()))[:60],
         'platforms': ';'.join(sorted(g.platform_id.dropna().unique())), 'char_keys': keys[:400]}
    for k, pat in FLAGS.items(): r[k] = int(bool(re.search(pat, text, re.I)))
    r['n_patients_guess'] = g.n_pat.iloc[0]
    rows.append(r)
F = pd.DataFrame(rows).merge(cat_df[['acc','title','type','date']], left_on='series', right_on='acc', how='left').drop(columns='acc')
F['clinical_score'] = F[['nephritis','sledai','response','proteinuria','histology','timepoint','treatment','autoantibody']].sum(axis=1)
F.sort_values(['clinical_score','n_samples'], ascending=False).to_csv(f'{O}/geo_series_label_flags.csv', index=False)
print(len(idx), 'samples;', len(F), 'series'); print(F.organism.value_counts().head())
