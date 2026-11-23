"""Audit r120: re-read Mackensen 2022 Nat Med source data (MOESM4 Fig 3D/3F, MOESM5 Fig 4B).
Usage: python3 -I r118_mackensen_recheck.py DIR_WITH_MOESM4_AND_MOESM5 (download: experiments/pc-reservoir/download.sh)"""
import sys, openpyxl
D = sys.argv[1]
ws = openpyxl.load_workbook(D + '/MOESM5.xlsx', data_only=True).worksheets[0]
rows = list(ws.iter_rows(values_only=True))
i0 = next(i for i, r in enumerate(rows) if r[0] == 'SLEDAI')
print('Fig 4B clinical anti-dsDNA (IU/mL), columns = patients 1-5')
curves = {p: [] for p in range(1, 6)}
for r in rows[i0 + 2:]:
    if r[7] is None: break
    for p in range(1, 6):
        v = r[7 + p]
        if v is not None: curves[p].append((r[7], v))
for p, c in curves.items():
    post = [x for x in c if x[0] > 0]
    q = len(post) >= 3 and max(m for m, _ in post) * 30 >= 90
    print(f'P{p}', c, 'post points', len(post), 'qualifies(>=3 post, last>=90d)', q,
          'month3/baseline', round(c[1][1] / c[0][1], 4) if len(c) > 1 else None)
