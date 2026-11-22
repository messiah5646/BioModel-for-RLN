# resp_table.csv: CRR responders per arm per visit (20 add-on LN RCTs)
Rebuild: `python3 -I resp_table/build_resp_table.py ../../data/real/ctgov/lupus_with_results.json.gz resp_table.csv`
then `python3 -I resp_table/add_publications.py` (adds 3 published week-104 rows -> resp_table_plus_pub.csv).
Check: `python3 -I resp_table/compare_q17.py` (reproduces result86.txt) and `resp_table/check_frozen.py` (round 90 + W1).
Both CSVs are committed with `git add -f` (*.csv is git-ignored). Round 107, PREREG-Q17.md, RESULT-Q17.txt.
