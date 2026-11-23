"""Round 118 harvest record. Queries (NCBI esearch db=pmc, open access filter), 2026-11-23:
 1 (rituximab OR obinutuzumab) AND lupus AND anti-dsDNA AND months   -> 2466 hits, top 200 by relevance
 2 ("CAR T" OR "CAR-T" OR "chimeric antigen receptor") AND lupus AND anti-dsDNA -> 462, top 200
 3 daratumumab AND lupus AND anti-dsDNA -> 145 (all)
 4 (blinatumomab OR teclistamab OR "bispecific") AND lupus AND anti-dsDNA -> 216, top 200
 Union 558 PMC ids (candidate_pmc_ids.txt). Every <table-wrap> scanned for (a) rows starting with an
 anti-dsDNA label carrying >=4 numeric cells (scan_rows.txt) and (b) anti-dsDNA columns with >=4 time-labelled
 numeric rows (scan_cols.txt). The scan code is the two inline scripts reproduced in RESULT-P6.md."""
