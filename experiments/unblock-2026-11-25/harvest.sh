#!/bin/bash
# Round 137 open-data harvest: every file goes through fetch.sh; failures are logged, not retried forever.
cd "$(dirname "$0")/../.."
L=experiments/unblock-2026-11-25/harvest_list.tsv; M=experiments/unblock-2026-11-25/MANIFEST-harvest.tsv
F=experiments/unblock-2026-11-25/harvest_failures.tsv; O=data/real/harvest-r137
while IFS=$'\t' read -r name url; do
  [ -s "$O/$name" ] && grep -q "	$name	" "$M" 2>/dev/null && continue
  bash scripts/provenance/fetch.sh "$url" "$O/$name" "$M" >/dev/null 2>&1 || printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$name" "$url" >> "$F"
done < "$L"
echo HARVEST_DONE
