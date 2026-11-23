#!/usr/bin/env bash
# P5 push 2: ILLUMINATE-1/2 preprocessed whole-blood expression + series matrices (open GEO). Usage: download.sh OUTDIR
set -euo pipefail
OUT=${1:?outdir}; mkdir -p "$OUT"; cd "$OUT"
B=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn
for g in GSE88885 GSE88886; do
  [ -s $g.matrix.gz ] || curl -s --retry 3 -o $g.matrix.gz $B/$g/matrix/${g}_series_matrix.txt.gz
done
[ -s GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz ] || curl -s --retry 3 --max-time 3600 -o GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz $B/GSE88885/suppl/GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz
[ -s GSE88886_ILLUMINATE2_pd_preprocessed.txt.gz ] || curl -s --retry 3 --max-time 3600 -o GSE88886_ILLUMINATE2_pd_preprocessed.txt.gz $B/GSE88886/suppl/GSE88886_ILLUMINATE2_pd_preprocessed.txt.gz
for f in *.gz; do python3 -I -c "import sys;h=open(sys.argv[1],'rb').read(2);print(sys.argv[1], 'OK' if h==b'\x1f\x8b' else 'NOT SERVED')" "$f"; done
sha256sum *.gz > SHA256SUMS
