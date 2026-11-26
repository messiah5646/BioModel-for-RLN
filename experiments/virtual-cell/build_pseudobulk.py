"""VC-1 (round 142) pseudobulk builder. Allowed before the 1-h mark: builds tables only, scores nothing.

Outputs (experiments/virtual-cell/work/, gitignored large files):
  units.csv            one row per (patient, compartment, cell_type, timepoint): n_cells, eligible (>= 30 cells)
  unit_sums.npy        full-gene raw count sums per unit (units x genes), float64
  genes.txt            gene ids (var_names) shared by all files
  hvg.txt              2,000 HVGs chosen on Base-timepoint eligible units only (PREREG section 4)
  cells_hvg.npz        per-cell sparse raw counts on the HVGs (for cell-bootstrap augmentation)
  cells.csv            per-cell unit index and library size (all-gene total counts)
"""
import sys, os, re, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

ROOT = os.path.dirname(os.path.abspath(__file__))
H = os.path.join(ROOT, '..', '..', 'data', 'real', 'harvest-r137')
OUT = os.path.join(ROOT, 'work'); os.makedirs(OUT, exist_ok=True)
FILES = {'T': 'CD3xCD19_Tcells.h5ad', 'B': 'CD3xCD19_Bcells.h5ad', 'Mono': 'CD3xCD19_Monocytes.h5ad',
         'Myeloid': 'CD3xCD19_Myeloid.h5ad', 'LDG': 'CD3xCD19_LDG.h5ad'}
MIN_CELLS = 30; NHVG = 2000; CHUNK = 20000
EXCL = re.compile(r'^(MT-|RPL|RPS|HB[ABDEGMQZ]\d?$|HBA|HBB|XIST$|TSIX$|RPS4Y|DDX3Y$|UTY$|KDM5D$|EIF1AY$|USP9Y$|ZFY$|NLGN4Y$|TXLNGY$)')


def raw_of(a):
    return a.raw.X if a.raw is not None else a.X


def main():
    units, sums, cell_unit, cell_lib, cell_blocks = [], [], [], [], []
    genes = None; uidx = {}
    for comp, fn in FILES.items():
        a = ad.read_h5ad(os.path.join(H, fn), backed='r')
        g = list(a.raw.var_names if a.raw is not None else a.var_names)
        if genes is None: genes = g
        assert g == genes, f'gene order differs in {fn}'
        o = a.obs[['donor_id', 'timepoint', 'cell_type', 'dose']].astype(str).reset_index(drop=True)
        key = o['donor_id'] + '|' + comp + '|' + o['cell_type'] + '|' + o['timepoint']
        for k in key.unique():
            if k not in uidx:
                p, c, ct, t = k.split('|')
                uidx[k] = len(units)
                units.append(dict(unit=len(units), patient=p, compartment=comp, cell_type=ct, timepoint=t,
                                  dose=o.loc[key.eq(k).idxmax(), 'dose'], n_cells=0))
                sums.append(np.zeros(len(genes)))
        ku = key.map(uidx).to_numpy()
        X = raw_of(a)
        for s in range(0, a.n_obs, CHUNK):
            e = min(s + CHUNK, a.n_obs)
            M = sp.csr_matrix(X[s:e])
            uu = ku[s:e]
            lib = np.asarray(M.sum(1)).ravel()
            for u in np.unique(uu):
                m = uu == u
                sums[u] += np.asarray(M[m].sum(0)).ravel()
                units[u]['n_cells'] += int(m.sum())
            cell_unit.append(uu); cell_lib.append(lib); cell_blocks.append((comp, s, e))
            print(comp, e, '/', a.n_obs, flush=True)
    U = pd.DataFrame(units); U['eligible'] = U.n_cells >= MIN_CELLS
    S = np.vstack(sums)
    U.to_csv(os.path.join(OUT, 'units.csv'), index=False)
    np.save(os.path.join(OUT, 'unit_sums.npy'), S)
    open(os.path.join(OUT, 'genes.txt'), 'w').write('\n'.join(genes))
    # HVGs on Base eligible units only (no outcome timepoints)
    base = U.index[(U.timepoint == 'Base') & U.eligible].to_numpy()
    L = np.log1p(S[base] / S[base].sum(1, keepdims=True) * 1e4)
    keep = np.array([not EXCL.match(x) for x in genes]) & (S[base].sum(0) > 0)
    var = np.where(keep, L.var(0), -1)
    hvg = np.argsort(-var)[:NHVG]
    open(os.path.join(OUT, 'hvg.txt'), 'w').write('\n'.join(str(i) for i in hvg))
    # per-cell HVG counts for bootstrap augmentation
    blocks = []
    for comp, fn in FILES.items():
        a = ad.read_h5ad(os.path.join(H, fn), backed='r'); X = raw_of(a)
        for s in range(0, a.n_obs, CHUNK):
            blocks.append(sp.csr_matrix(X[s:min(s + CHUNK, a.n_obs)])[:, hvg].astype(np.float32))
    sp.save_npz(os.path.join(OUT, 'cells_hvg.npz'), sp.vstack(blocks).tocsr())
    pd.DataFrame({'unit': np.concatenate(cell_unit), 'lib': np.concatenate(cell_lib)}).to_csv(
        os.path.join(OUT, 'cells.csv'), index=False)
    print(json.dumps({'units': len(U), 'eligible': int(U.eligible.sum()), 'cells': int(U.n_cells.sum()),
                      'base_units_for_hvg': int(len(base))}))


if __name__ == '__main__':
    main()
