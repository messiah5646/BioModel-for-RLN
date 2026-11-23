# PREREG: audit round 111 (rounds 91-110), committed before any re-run

Scope (rule 10, queue item P1).
A. Reproduction (scripts run in a scratch copy of experiments/lumen-bench, committed files untouched):
  A1 build_resp_table.py on data/real/ctgov/lupus_with_results.json.gz -> diff vs committed resp_table.csv. Pass: identical rows; counts equal.
  A2 compare_q17.py (P1-P3 vs result86.txt). Pass: prints PASS (tolerances as in PREREG-Q17: lor/se 0.01, IV-MSE 0.002).
  A3 check_frozen.py (round 90 + W1 frozen forecasts). Pass: reproduces round-90 MSE 0.055 vs carry-forward 0.242 within 0.005.
B. Claim ledger: every Supported/Established claim in rounds 91-110 with n, pre-registered (Y/N), test count.
  Downgrade rule: Supported/Established on n < 10 patients/trials, or not pre-registered, -> Emerging.
  Additionally flag any claim whose load-bearing number cannot be traced to committed code/output.
C. Verdict: audit PASSES if A1-A3 pass and no claim is found mislabelled after B; otherwise FAILS with the list.
