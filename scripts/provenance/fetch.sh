#!/bin/bash
# Rule 12b helper: download one file and append a provenance line to a committed manifest.
# Usage: fetch.sh URL OUTFILE MANIFEST.tsv
# Manifest columns: fetch_utc  url  outfile  bytes  md5  http_last_modified
set -euo pipefail
URL=$1; OUT=$2; MAN=$3
HDR=$(mktemp)
# NCBI returns transient 403s (seen twice in audit r125); retry up to 4 times, 30 s apart
for i in 1 2 3 4; do curl -sSfL -D "$HDR" -o "$OUT" "$URL" && break; [ $i -eq 4 ] && exit 22; sleep 30; done
LM=$(grep -i '^last-modified:' "$HDR" | tail -1 | cut -d' ' -f2- | tr -d '\r' || true)
rm -f "$HDR"
[ -f "$MAN" ] || printf 'fetch_utc\turl\toutfile\tbytes\tmd5\thttp_last_modified\n' > "$MAN"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$URL" "$(basename "$OUT")" \
  "$(stat -c %s "$OUT")" "$(md5sum "$OUT" | cut -d' ' -f1)" "${LM:-NA}" >> "$MAN"
tail -1 "$MAN"
