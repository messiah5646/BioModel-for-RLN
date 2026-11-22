# Open-data harvest (Q0, 2026-11-22)

Policy: open data only (no account, no data-use agreement). See AUTONOMOUS-QUEUE.md.

| File (data/real/catalogues/) | What |
|---|---|
| geo_lupus_series.csv | 597 human/mouse lupus GEO series (2026-11-21 search) |
| geo_sample_index.tsv.gz | 62,446 GEO samples with all sample characteristics (from series-matrix headers) |
| geo_series_label_flags.csv | per-series flags: nephritis, SLEDAI/activity, response, proteinuria, histology, timepoint, treatment, autoantibody; clinical_score = sum |
| geo_clinical_files_manifest.csv | every matrix/suppl file (with size) of the 28 series with clinical_score >= 3 (152 GB total, 0.49 GB of it series matrices) |
| arrayexpress_lupus_nongeo.csv | ArrayExpress lupus entries that are not GEO mirrors (keyword hits include unrelated canine studies) |
| cellxgene_lupus_datasets.json | all 6 CELLxGENE datasets tagged lupus / lupus nephritis, with H5AD URLs |

Re-create in a fresh workspace (about 10 min, ~0.5 GB):

    python3 -I scripts/harvest/geo_headers.py data/real/catalogues/geo_lupus_series.csv /home/claude/data/geo_headers
    python3 -I scripts/harvest/build_index.py /home/claude/data/geo_headers data/real/catalogues/geo_lupus_series.csv data/real/catalogues
    # expression matrices of the clinically labelled series:
    python3 -c "import pandas as pd;m=pd.read_csv('data/real/catalogues/geo_clinical_files_manifest.csv');m=m[m.dir=='matrix'];print('\n'.join(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}/matrix/{f}' for a,f in zip(m.series,m.file)))" | xargs -P4 -n1 curl -s -O

Not downloaded on purpose: raw CEL/FASTQ supplements (~150 GB, mostly ILLUMINATE GSE88884-7). Expression for RNA-seq series lives in suppl/ and is fetched per test.
Not open (checked 2026-11-22): every AMP RA/SLE file in the ARK Portal (Synapse) carries the controlled-access requirement; ImmPort SDY997 and dbGaP need accounts.
