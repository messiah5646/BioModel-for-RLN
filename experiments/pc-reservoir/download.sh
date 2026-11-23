#!/bin/bash
# P6 data acquisition (round 118). Open sources only.
set -e; D=${1:-/home/claude/data/p6}; mkdir -p $D/mackensen2022 $D/feng2025
# Mackensen 2022 Nat Med (doi 10.1038/s41591-022-02017-5) source data
for n in 2 3 4 5 6 7 8; do curl -s -o $D/mackensen2022/MOESM$n.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-022-02017-5/MediaObjects/41591_2022_2017_MOESM${n}_ESM.xlsx"; done
curl -s -o $D/mackensen2022/MOESM1.pdf "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-022-02017-5/MediaObjects/41591_2022_2017_MOESM1_ESM.pdf"
# Feng 2025 Nat Med (doi 10.1038/s41591-025-03937-8) supplementary tables
curl -s -o $D/feng2025/MOESM3.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-025-03937-8/MediaObjects/41591_2025_3937_MOESM3_ESM.xlsx"
# PMC open-access table scan: NCBI esearch (4 queries, see scan_pmc.py) -> candidate_pmc_ids.txt (558), efetch XML
