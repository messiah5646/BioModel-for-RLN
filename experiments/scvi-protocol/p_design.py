"""Q36 design inspection (obs only, no expression). Usage: python3 -I p_design.py <h5ad> <out_dir>"""
import sys, h5py, numpy as np, pandas as pd
H, O = sys.argv[1:3]
f = h5py.File(H, 'r'); g = f['obs']
def col(k):
    x = g[k]
    if isinstance(x, h5py.Group): return pd.Categorical.from_codes(x['codes'][:], x['categories'][:].astype(str))
    return x[:].astype(str) if x.dtype.kind in 'SO' else x[:]
print('obs keys', list(g.keys())); print('raw keys', list(f['raw'].keys()) if 'raw' in f else None, 'X', type(f['X']))
cols = ['Processing_Cohort', 'disease', 'disease_state', 'donor_id', 'library_uuid', 'sample_uuid', 'author_cell_type']
M = pd.DataFrame({k: col(k) for k in cols}); M.to_pickle(f'{O}/obs.pkl')
print(M.Processing_Cohort.value_counts())
for c in sorted(M.Processing_Cohort.unique()):
    m = M[M.Processing_Cohort == c]
    print('cohort', c, m.groupby('disease_state').donor_id.nunique().to_dict(), 'libs', m.library_uuid.nunique(), 'cells', len(m))
c3 = M[M.disease_state.isin(['flare', 'treated'])]
st = c3.groupby('donor_id').disease_state.agg(lambda s: tuple(sorted(set(s))))
print('flare/treated donors', st.value_counts().to_dict(), 'cohorts', c3.Processing_Cohort.unique().tolist())
both = st[st.apply(len) == 2].index
pr = c3[c3.donor_id.isin(both)]
print('paired donors', len(both)); print(pr.groupby(['donor_id', 'disease_state']).size().unstack())
print('paired: libraries shared by flare & treated of same donor', sum(len(set(pr[(pr.donor_id==d)&(pr.disease_state=='flare')].library_uuid) & set(pr[(pr.donor_id==d)&(pr.disease_state=='treated')].library_uuid)) > 0 for d in both))
print(M.author_cell_type.value_counts().to_dict())
br = M.groupby('donor_id').Processing_Cohort.nunique(); print('bridging donors', int((br > 1).sum()))
