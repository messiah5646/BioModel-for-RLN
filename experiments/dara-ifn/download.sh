#!/usr/bin/env bash
# P2b push 2: download GSE294789 gene-expression h5 + TCR/BCR contigs (open GEO). Usage: download.sh OUTDIR
set -euo pipefail
OUT=${1:?outdir}; mkdir -p "$OUT"; cd "$OUT"
curl -s --max-time 120 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294789&targ=all&form=text&view=brief" | tr -d "\r" > soft.txt
awk '/^\^SAMPLE/{s=$3} /^!Sample_title/{t=substr($0,index($0,"=")+2)} /^!Sample_supplementary_file_1/{f=$3; sub(/^ftp:/,"https:",f); print s"\t"t"\t"f}' soft.txt > manifest.tsv
cut -f3 manifest.tsv | while read -r u; do
  f=$(basename "$u"); [ -s "$f" ] || curl -s --retry 3 --max-time 1800 -o "$f" "$u"
done
sha256sum *.h5 *.csv.gz > SHA256SUMS
wc -l manifest.tsv SHA256SUMS
# report files GEO did not serve (403 pages saved instead of data)
python3 -I -c "
import glob
for f in sorted(glob.glob('*.gz')+glob.glob('*.h5')):
    h=open(f,'rb').read(4)
    if not (h[:2]==b'\x1f\x8b' if f.endswith('.gz') else h==b'\x89HDF'): print('NOT SERVED', f)
" | tee NOT_SERVED.txt
