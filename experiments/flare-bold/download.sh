#!/bin/bash
# Round 122 data: GSE92776 series matrix + platform annotation (GEO, open)
set -e; D=${1:-/home/claude/data/geo/GSE92776}; mkdir -p $D; cd $D
curl -sSfO https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92776/matrix/GSE92776_series_matrix.txt.gz
curl -sSf "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL22843&targ=self&form=text&view=data" -o GPL22843.txt
md5sum GSE92776_series_matrix.txt.gz GPL22843.txt
