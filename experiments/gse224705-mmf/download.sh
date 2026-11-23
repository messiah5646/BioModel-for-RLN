#!/bin/bash
# Usage: bash download.sh <data_dir>. Rebuilds the inputs run.py / q25_trajectory.py read (GSE224705, open GEO).
set -e; D=$1; mkdir -p $D
curl -sS -o $D/sm.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224705/matrix/GSE224705_series_matrix.txt.gz
curl -sS -o $D/gpl.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz
python3 -I - $D <<'PY'
import sys, gzip, csv
D = sys.argv[1]; rows = {}; gsm = None
for l in gzip.open(f'{D}/sm.txt.gz', 'rt'):
    if l.startswith('!Sample_geo_accession'): gsm = [x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
    elif l.startswith('!Sample_characteristics_ch1'):
        v = [x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
        k = v[0].split(': ', 1)[0]; rows[k] = [x.split(': ', 1)[1] if ': ' in x else '' for x in v]
    elif l.startswith('!series_matrix_table_begin'): break
keys = list(rows)
with open(f'{D}/chars.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['gsm'] + keys)
    for i, g in enumerate(gsm): w.writerow([g] + [rows[k][i] for k in keys])
print('samples', len(gsm), 'fields', keys)
PY
