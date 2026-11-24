"""Audit r135: blind flag scan of the abstracts committed in newtrials/raw/pubmed_efetch.xml (836f753).
Prints FLAGS ONLY (no abstract text, no numbers) so the auditor stays blind to outcomes."""
import csv, re, sys, os, xml.etree.ElementTree as ET
N = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lumen-bench', 'newtrials')
root = ET.parse(os.path.join(N, 'raw', 'pubmed_efetch.xml')).getroot()
ab = {a.findtext('.//MedlineCitation/PMID'): ' '.join(''.join(x.itertext()) for x in a.findall('.//Abstract/AbstractText'))
      for a in root.findall('.//PubmedArticle')}
EARLY = re.compile(r'\b(week|wk)s?\s*(12|16|20|24|26)\b|\b(12|16|20|24|26)[- ]weeks?\b|\bmonth\s*(3|6)\b|\b(3|6)[- ]months?\b', re.I)
NUM = re.compile(r'\d+(\.\d+)?\s*%|\b\d+\s*/\s*\d+\b|\bn\s*=\s*\d+', re.I)
w = csv.writer(sys.stdout, delimiter='\t'); w.writerow(['cid', 'indication', 'pmid', 'abstract_committed', 'numeric', 'early_numeric_sentence'])
for r in csv.DictReader(open(os.path.join(N, 'candidates_frozen.csv'))):
    t = ab.get(r['pmid'], '') if r['pmid'] else ''
    es = any(EARLY.search(s) and NUM.search(s) for s in re.split(r'(?<=[.;])\s', t))
    w.writerow([r['cid'], r['indication'], r['pmid'] or '-', int(bool(t)), int(bool(NUM.search(t))), int(es)])
