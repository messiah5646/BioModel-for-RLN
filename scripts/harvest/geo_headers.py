"""Stream only the metadata header of every GEO series matrix (stops at the expression table).
Usage: python3 -I geo_headers.py <catalogue.csv> <out_dir>. Writes one TSV of sample characteristics per series."""
import sys, os, re, gzip, zlib, urllib.request, csv, concurrent.futures as cf
cat, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
accs = [r['acc'] for r in csv.DictReader(open(cat))]
def base(a): return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}/matrix/"
def head(url):
    d = zlib.decompressobj(16 + zlib.MAX_WBITS); buf = b''
    with urllib.request.urlopen(url, timeout=120) as r:
        while True:
            c = r.read(1 << 16)
            if not c: break
            buf += d.decompress(c)
            if b'!series_matrix_table_begin' in buf: break
    return buf.split(b'!series_matrix_table_begin')[0].decode('utf-8', 'replace')
def one(a):
    dst = f"{out}/{a}.tsv"
    if os.path.exists(dst): return a, 'cached'
    try:
        names = re.findall(r'href="(GSE[^"]+_series_matrix\.txt\.gz)"', urllib.request.urlopen(base(a), timeout=60).read().decode())
        rows = []
        for n in sorted(set(names)):
            h = head(base(a) + n); lines = [l.split('\t') for l in h.splitlines() if l.startswith('!Sample_')]
            S = {}
            for l in lines:
                k = l[0][8:]; v = [x.strip('"') for x in l[1:]]
                if k == 'characteristics_ch1':
                    for j, x in enumerate(v):
                        if ': ' in x:
                            kk, vv = x.split(': ', 1); S.setdefault('ch:' + kk.strip().lower(), [''] * len(v))[j] = vv
                elif k in ('geo_accession', 'title', 'source_name_ch1', 'organism_ch1', 'platform_id', 'molecule_ch1', 'library_strategy'):
                    S.setdefault(k, v)
            n_s = len(S.get('geo_accession', []))
            for j in range(n_s): rows.append({'series': a, 'matrix': n, **{k: (vals[j] if j < len(vals) else '') for k, vals in S.items()}})
        keys = sorted({k for r in rows for k in r}, key=lambda k: (k.startswith('ch:'), k))
        with open(dst, 'w', newline='') as f:
            w = csv.DictWriter(f, keys, delimiter='\t'); w.writeheader(); w.writerows(rows)
        return a, len(rows)
    except Exception as e:
        return a, f'ERR {type(e).__name__}: {e}'[:120]
with cf.ThreadPoolExecutor(6) as ex:
    for a, r in ex.map(one, accs): print(a, r, flush=True)
