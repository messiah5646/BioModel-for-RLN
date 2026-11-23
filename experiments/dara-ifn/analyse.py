"""P2b analysis (push 3) per experiments/dara-ifn/PREREG.md. Usage: python analyse.py DATADIR OUTDIR
DATADIR = output of download.sh (h5 + TCR/BCR contig csv.gz + manifest.tsv)."""
import sys, os, glob, json, gzip, itertools
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy.stats import spearmanr

D, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
G = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']
QC = ['MT-']  # mito prefix

# Library -> patient-samples, hand-coded from GEO titles (manifest.tsv). (patient, week)
LIB = {
 'GSM8919721': [(1, 0)], 'GSM8919724': [(2, 0)], 'GSM8919727': [(1, 9)], 'GSM8919730': [(2, 9)],
 'GSM8919733': [(3, 0)], 'GSM8919736': [(3, 9)], 'GSM8919739': [(4, 0)], 'GSM8919742': [(5, 0), (7, 0)],
 'GSM8919745': [(4, 9)], 'GSM8919748': [(1, 36), (2, 36)], 'GSM8919751': [(5, 9), (7, 9)],
 'GSM8919754': [(3, 36)], 'GSM8919757': [(9, 0)], 'GSM8919760': [(5, 36), (10, 0)],
 'GSM8919763': [(7, 36), (9, 9)], 'GSM8919766': [(9, 9)], 'GSM8919769': [(10, 9)], 'GSM8919772': [(9, 36)],
 'GSM8919775': [(8, 36), (10, 36)],
}
# anti-dsDNA IU/ml, source data Fig 2a1 (paper SLE #N == GEO Patient N, assumption A1)
DS = {1: (331, 73.4, 94.6), 2: (1769, 201, 127.2), 3: (64.9, 27.8, 31.3), 4: (167.1, 61.7, 66.2), 5: (255, 89.7, 107.3),
      6: (87.3, 22, 27.8), 7: (59.2, 22.3, 25), 8: (165.5, 54.4, 20.4), 9: (144.7, 114, 197.4), 10: (657, 128.6, 392)}  # w0, w8, w36

man = pd.read_csv(os.path.join(D, 'manifest.tsv'), sep='\t', header=None, names=['gsm', 'title', 'url'])
man['file'] = man.url.map(os.path.basename)
def vdj_file(gsm_expr, kind):
    i = man.index[man.gsm == gsm_expr][0]
    for j in (i + 1, i + 2):
        if j in man.index and kind in man.file[j]: return os.path.join(D, man.file[j])
    return None  # not deposited (GSM8919775 has no TCR)

def load_vdj(path, chains):
    if path is None or open(path, 'rb').read(2) != b'\x1f\x8b':
        return pd.DataFrame(columns=['barcode', 'chain', 'cdr3_nt'])
    v = pd.read_csv(path)
    v = v[(v.productive.astype(str).str.lower() == 'true') & (v.high_confidence.astype(str).str.lower() == 'true') & v.chain.isin(chains)]
    return v

def read_h5(path):
    m = h5py.File(path, 'r')['matrix']
    X = sp.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=m['shape'][:])  # features x cells
    name = np.array(m['features']['name']).astype(str); ft = np.array(m['features']['feature_type']).astype(str)
    bc = np.array(m['barcodes']).astype(str)
    return X, name, ft, bc

# ---- per library: cell table ----
cells = []; skipped = []
for gsm in LIB:
    h5 = glob.glob(os.path.join(D, gsm + '_*.h5'))[0]
    if open(h5, 'rb').read(4) != b'\x89HDF':
        skipped.append(gsm); continue  # GEO returned 403 / not HDF5
    X, name, ft, bc = read_h5(h5)
    ge = ft == 'Gene Expression'
    Xg = X[ge].tocsr(); gname = name[ge]
    tot = np.asarray(Xg.sum(0)).ravel(); ngen = np.asarray((Xg > 0).sum(0)).ravel()
    mt = np.asarray(Xg[np.char.startswith(gname, 'MT-')].sum(0)).ravel() / np.maximum(tot, 1)
    sc = np.mean([np.log1p(Xg[np.where(gname == g)[0][0]].toarray().ravel() / np.maximum(tot, 1) * 1e4) for g in G], axis=0)
    ab = X[ft == 'Antibody Capture'].tocsr(); abn = name[ft == 'Antibody Capture']
    h1 = ab[np.where(abn == 'Hashtag1')[0][0]].toarray().ravel() if 'Hashtag1' in abn else np.zeros(len(bc))
    h2 = ab[np.where(abn == 'Hashtag2')[0][0]].toarray().ravel() if 'Hashtag2' in abn else np.zeros(len(bc))
    t = load_vdj(vdj_file(gsm, 'TCR'), ['TRB']); b = load_vdj(vdj_file(gsm, 'BCR'), ['IGH'])
    trb = t.groupby('barcode').cdr3_nt.apply(lambda s: '|'.join(sorted(set(s))))
    df = pd.DataFrame(dict(lib=gsm, bc=bc, tot=tot, ngen=ngen, mt=mt, ifn6=sc, h1=h1, h2=h2))
    df['trb'] = df.bc.map(trb); df['has_t'] = df.trb.notna(); df['has_b'] = df.bc.isin(set(b.barcode))
    df['ctype'] = np.where(df.has_t & ~df.has_b, 'T', np.where(df.has_b & ~df.has_t, 'B', 'other'))
    df['no_tcr_file'] = len(t) == 0
    if len(t) == 0:  # NOT PRE-REGISTERED: CD3E/CD3D > 0 and no BCR -> 'T_rna'
        cd3 = (Xg[np.where(gname == 'CD3E')[0][0]].toarray().ravel() + Xg[np.where(gname == 'CD3D')[0][0]].toarray().ravel()) > 0
        df['ctype'] = np.where(cd3 & ~df.has_b, 'T_rna', df.ctype)
    igh = b.groupby('barcode').cdr3_nt.apply(lambda s: '|'.join(sorted(set(s)))) if len(b) else pd.Series(dtype=str)
    df['igh'] = df.bc.map(igh)
    df = df[df.tot > 500]
    cells.append(df)
cells = pd.concat(cells, ignore_index=True)
cells['patient'] = np.nan; cells['week'] = np.nan
single = {g: v[0] for g, v in LIB.items() if len(v) == 1}
for g, (p, w) in single.items():
    cells.loc[cells.lib == g, ['patient', 'week']] = (p, w)

# ---- A2: TCR clonotype demultiplexing of pooled libraries ----
def clones(mask):
    s = set()
    for k in cells.loc[mask & (cells.ctype == 'T'), 'trb'].dropna(): s.update(k.split('|'))
    return s
ref = {p: clones(cells.patient == p) for p in range(1, 11)}
pooled = {g: v for g, v in LIB.items() if len(v) == 2}
# patients with no single-patient library get a reference bootstrapped from clones shared between two of their pools
def bootstrap(p, libs_with_p, known_partner):
    s = None
    for g, partner in zip(libs_with_p, known_partner):
        m = (cells.lib == g) & (cells.ctype == 'T')
        c = set()
        for k in cells.loc[m, 'trb'].dropna():
            parts = k.split('|')
            if partner is None or not any(x in ref[partner] for x in parts): c.update(parts)
        s = c if s is None else s & c
    return s
# P5: pools 742 (with P7), 751 (with P7), 760 (with P10, known). P7: 742, 751, 763 (with P9, known).
s5 = bootstrap(5, ['GSM8919760', 'GSM8919742'], [10, None]) | bootstrap(5, ['GSM8919760', 'GSM8919751'], [10, None])
s7 = bootstrap(7, ['GSM8919763', 'GSM8919742'], [9, None]) | bootstrap(7, ['GSM8919763', 'GSM8919751'], [9, None])
ref[5] = s5 - s7; ref[7] = s7 - s5
demux = []
for g, ((pa, wa), (pb, wb)) in pooled.items():
    m = (cells.lib == g) & (cells.ctype == 'T')
    for i in cells.index[m]:
        k = cells.at[i, 'trb']; parts = k.split('|')
        ia = any(x in ref[pa] for x in parts); ib = any(x in ref[pb] for x in parts)
        if ia and not ib: cells.loc[i, ['patient', 'week']] = (pa, wa)
        elif ib and not ia: cells.loc[i, ['patient', 'week']] = (pb, wb)
    # hashtag cross-check: which hashtag dominates each clonotype-assigned patient
    sub = cells[(cells.lib == g) & cells.patient.notna() & (cells.ctype == 'T')].copy()
    sub['ht'] = np.where(sub.h1 > sub.h2, 1, np.where(sub.h2 > sub.h1, 2, 0))
    tab = sub.groupby(['patient', 'ht']).size().unstack(fill_value=0)
    maj = {int(p): int(r.drop(0, errors='ignore').idxmax()) for p, r in tab.iterrows() if r.drop(0, errors='ignore').sum() > 0}
    conc = np.mean([sub.ht[sub.patient == p].eq(h).mean() for p, h in maj.items()]) if maj else np.nan
    # assign B cells of pooled libraries by the learned hashtag mapping (only if mapping is one-to-one)
    if len(maj) == 2 and len(set(maj.values())) == 2:
        inv = {h: p for p, h in maj.items()}; wk = {pa: wa, pb: wb}
        mb = (cells.lib == g) & (cells.ctype == 'B')
        htb = np.where(cells.loc[mb, 'h1'] > cells.loc[mb, 'h2'], 1, np.where(cells.loc[mb, 'h2'] > cells.loc[mb, 'h1'], 2, 0))
        for i, h in zip(cells.index[mb], htb):
            if h in inv: cells.loc[i, ['patient', 'week']] = (inv[h], wk[inv[h]])
    demux.append(dict(lib=g, pair=f'P{pa}w{wa}+P{pb}w{wb}', t_cells=int(m.sum()),
                      assigned_a=int(((cells.lib == g) & (cells.patient == pa) & (cells.ctype == 'T')).sum()),
                      assigned_b=int(((cells.lib == g) & (cells.patient == pb) & (cells.ctype == 'T')).sum()),
                      hashtag_map=json.dumps(maj), hashtag_concordance=round(float(conc), 3)))
demux = pd.DataFrame(demux)

# ---- QC-residualised IFN6 (per-cell linear model on all T cells) ----
T = cells[cells.ctype == 'T']
Z = np.column_stack([np.ones(len(T)), np.log(T.tot), np.log(T.ngen), T.mt])
beta, *_ = np.linalg.lstsq(Z, T.ifn6.values, rcond=None)
cells['ifn6_res'] = np.nan
cells.loc[T.index, 'ifn6_res'] = T.ifn6.values - Z @ beta + T.ifn6.mean()

# ---- per patient-sample summary ----
a = cells[cells.patient.notna()]
ps = a.groupby(['patient', 'week']).apply(lambda g: pd.Series(dict(
    nT=int((g.ctype == 'T').sum()), nB=int((g.ctype == 'B').sum()),
    ifn6=g.ifn6[g.ctype == 'T'].mean(), ifn6_res=g.ifn6_res[g.ctype == 'T'].mean()))).reset_index()
ps['bfrac'] = ps.nB / (ps.nB + ps.nT)
ps = ps[ps.nT >= 50]  # A2 drop rule
ps.to_csv(os.path.join(OUT, 'per_sample.csv'), index=False)

def ratios(col, wk):
    rows = []
    for p in sorted(ps.patient.unique()):
        b0 = ps[(ps.patient == p) & (ps.week == 0)]; bt = ps[(ps.patient == p) & (ps.week == wk)]
        if len(b0) and len(bt):
            ds = DS[int(p)][2 if wk == 36 else 1] / DS[int(p)][0]
            rows.append(dict(patient=int(p), m=float(bt[col].iloc[0] / b0[col].iloc[0]), ds_ratio=ds,
                             bfrac_ratio=float(bt.bfrac.iloc[0] / b0.bfrac.iloc[0]) if b0.bfrac.iloc[0] > 0 else np.nan))
    return pd.DataFrame(rows)

def exact_perm_p(x, y):
    r0 = spearmanr(x, y)[0]; n = len(x)
    rs = [spearmanr(x, np.array(y)[list(perm)])[0] for perm in itertools.permutations(range(n))]
    return float(np.mean(np.array(rs) >= r0 - 1e-12))

res = {}
for col in ['ifn6', 'ifn6_res']:
    r36 = ratios(col, 36); r9 = ratios(col, 9)
    rho36 = float(spearmanr(np.log(r36.ds_ratio), np.log(r36.m))[0]) if len(r36) >= 3 else np.nan
    rho9 = float(spearmanr(np.log(r9.ds_ratio), np.log(r9.m))[0]) if len(r9) >= 3 else np.nan
    n = len(r36)
    R1 = 'H1 supported' if (rho36 >= 0.5 and n >= 5) else ('H1 fails' if rho36 <= 0 else 'inconclusive')
    medm = float(r36.m.median()); medb = float(r36.bfrac_ratio.median())
    R2 = ('B-cell explanation fails' if (medm <= 0.80 and medb >= 0.8) else
          ('inconclusive (B fraction fell)' if medb < 0.8 else 'B-cell explanation not rejected'))
    res[col] = dict(n36=n, rho36=rho36, perm_p36=exact_perm_p(np.log(r36.ds_ratio).values, np.log(r36.m).values) if n >= 3 else None,
                    R1=R1, median_m36=medm, median_bfrac_ratio36=medb, R2=R2,
                    R3_frac_fall36=float((r36.m < 1).mean()), n9=len(r9), rho9=rho9, median_m9=float(r9.m.median()),
                    per_patient36=r36.round(3).to_dict('records'), per_patient9=r9.round(3).to_dict('records'))
res['skipped_libraries'] = skipped
res['qc_flip'] = res['ifn6']['R1'] != res['ifn6_res']['R1']
res['final_R1'] = 'inconclusive (QC flip)' if res['qc_flip'] else res['ifn6']['R1']
demux.to_csv(os.path.join(OUT, 'demux.csv'), index=False)
json.dump(res, open(os.path.join(OUT, 'result.json'), 'w'), indent=1)
print(demux.to_string()); print(ps.round(4).to_string()); print(json.dumps({k: v for k, v in res.items()}, indent=1))

# ================= NOT PRE-REGISTERED sensitivity: P10 week 36 via BCR clonotype -> hashtag ===============
# GSM8919775 (P8w36 + P10w36) has no TCR deposit, so A2 drops P10 from the primary analysis.
g = 'GSM8919775'
bref10 = set(); [bref10.update(k.split('|')) for k in cells.loc[(cells.patient == 10) & (cells.ctype == 'B'), 'igh'].dropna()]
L = cells[cells.lib == g].copy()
L['ht'] = np.where(L.h1 > L.h2, 1, np.where(L.h2 > L.h1, 2, 0))
Bm = L[(L.ctype == 'B') & L.igh.notna()]
hit = Bm.igh.map(lambda k: any(x in bref10 for x in k.split('|')))
sens = dict(p10_bcr_matched_B=int(hit.sum()))
if hit.sum() >= 5:
    vc = Bm[hit].ht.value_counts().drop(0, errors='ignore')
    h10 = int(vc.idxmax()) if len(vc) else 0
    sens['p10_hashtag'] = h10; sens['hashtag_purity_matched_B'] = float((Bm[hit].ht == h10).mean())
    Tp10 = L[(L.ctype == 'T_rna') & (L.ht == h10)]
    b0 = ps[(ps.patient == 10) & (ps.week == 0)]
    sens['p10_w36_nT'] = len(Tp10)
    if len(Tp10) >= 50 and len(b0):
        sens['p10_m36'] = float(Tp10.ifn6.mean() / b0.ifn6.iloc[0])
        r = ratios('ifn6', 36)
        r = pd.concat([r, pd.DataFrame([dict(patient=10, m=sens['p10_m36'], ds_ratio=DS[10][2] / DS[10][0])])])
        sens['rho36_with_p10'] = float(spearmanr(np.log(r.ds_ratio), np.log(r.m))[0]); sens['n_with_p10'] = len(r)
        sens['note'] = 'T cells here are CD3E/CD3D+ BCR- (not TCR-defined), so P10 w36 is not comparable in definition to w0'
res['NOT_PREREG_sensitivity_p10'] = sens
json.dump(res, open(os.path.join(OUT, 'result.json'), 'w'), indent=1)
print('SENSITIVITY', json.dumps(sens, indent=1))

# ================= NOT PRE-REGISTERED: label audit + reliable subset ===============
# TCR clonotype overlap between libraries (fraction of the smaller library's TRB CDR3nt set) as an identity check
# on the GEO sample titles.
tr = {}
for gsm in LIB:
    v = load_vdj(vdj_file(gsm, 'TCR'), ['TRB'])
    if len(v): tr[gsm] = set(v.cdr3_nt)
v766 = load_vdj(os.path.join(D, man.file[man.index[man.gsm == 'GSM8919768'][0]]), ['TRB']); tr['GSM8919766'] = set(v766.cdr3_nt)
labs = {g: '+'.join(f'P{p}w{w}' for p, w in LIB[g]) for g in tr}
ov = pd.DataFrame({labs[a]: {labs[b]: round(100 * len(tr[a] & tr[b]) / min(len(tr[a]), len(tr[b])), 1) for b in tr} for a in tr})
ov.to_csv(os.path.join(OUT, 'label_audit_tcr_overlap_pct.csv'))
# reliable: patients whose w0 and w36 sit in libraries confirmed by clonotype overlap (P1, P3, P9) + P2 via hashtag
L = cells[(cells.lib == 'GSM8919748') & (cells.ctype == 'T')]
p2 = L[L.h2 > L.h1]  # hashtag 2 = P2 (clonotype-learned map, concordance in demux.csv)
rows = []
for p in [1, 2, 3, 9]:
    b0 = ps[(ps.patient == p) & (ps.week == 0)]
    if p == 2: w36 = p2.ifn6.mean(); w36r = p2.ifn6_res.mean(); n = len(p2)
    else:
        bt = ps[(ps.patient == p) & (ps.week == 36)]; w36, w36r, n = bt.ifn6.iloc[0], bt.ifn6_res.iloc[0], int(bt.nT.iloc[0])
    rows.append(dict(patient=p, nT36=n, m36=w36 / b0.ifn6.iloc[0], m36_res=w36r / b0.ifn6_res.iloc[0], ds_ratio=DS[p][2] / DS[p][0]))
rel = pd.DataFrame(rows)
res['NOT_PREREG_reliable_subset'] = dict(table=rel.round(3).to_dict('records'),
    rho=float(spearmanr(np.log(rel.ds_ratio), np.log(rel.m36))[0]), rho_res=float(spearmanr(np.log(rel.ds_ratio), np.log(rel.m36_res))[0]),
    note='GEO titles of GSM8919760 and GSM8919763 contradict TCR overlap; P5/P7/P10 w36 assignments unreliable')
json.dump(res, open(os.path.join(OUT, 'result.json'), 'w'), indent=1)
print('RELIABLE', json.dumps(res['NOT_PREREG_reliable_subset'], indent=1))
