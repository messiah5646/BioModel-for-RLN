#!/usr/bin/env bash
# Q22 push 2: download GSE304878 scRNA matrices (open GEO). Usage: download.sh OUTDIR
set -euo pipefail
OUT=${1:?outdir}; mkdir -p "$OUT"; cd "$OUT"
curl -s --max-time 120 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304878&targ=all&form=text&view=brief" | tr -d "\r" > soft.txt
grep -o 'ftp://[^ ]*scRNA[^ ]*\|ftp://[^ ]*_\(PRE\|POST\)[0-9]_\(barcodes\|features\|matrix\)[^ ]*' soft.txt | sort -u | sed 's/^ftp:/https:/' > urls.txt
for pass in 1 2 3 4; do while read -r u; do f=$(basename "$u"); if [ -s "$f" ] && [ "$(head -c2 "$f" | od -An -tx1 | tr -d " ")" = "1f8b" ]; then continue; fi; sleep 3; curl -s --retry 3 --max-time 1800 -o "$f" "$u"; done < urls.txt; done
sha256sum *.gz > SHA256SUMS; wc -l urls.txt SHA256SUMS
for f in *.gz; do [ "$(head -c2 "$f" | od -An -tx1 | tr -d " ")" = "1f8b" ] || echo "NOT SERVED $f"; done
