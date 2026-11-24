"""RT-20261124-1 run 1: title-level screen of PubMed RCT records (titles, dates, registry IDs only;
abstracts are not parsed or printed). Removes records whose NCT is among the 60 tabulated trials
and titles that are clearly secondary (post hoc, pooled, extension, subgroup, PK, QoL, cost)."""
import re, csv, os, sys, xml.etree.ElementTree as ET
D = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(D, '..')
known = set()
for f in ['resp_table.csv', 'sle/sri4_table.csv', 'sle/bicla_table.csv']:
    for r in csv.DictReader(open(os.path.join(B, f))): known.add(r['nct'])
SEC = re.compile(r'post[- ]?hoc|pooled|extension|long-term|subgroup|sub-?analys|secondary analys|pharmacokinet|'
                 r'quality of life|cost|protocol|design and rationale|rationale and design|exploratory|biomarker|'
                 r'patient-reported|fatigue|vaccin|pregnan|meta-analys|retrospective|open-label extension|follow-up|'
                 r'transcript|gene expression|interferon signature|cognit|osteopor|bone|lipid|atheroscl|cardiovasc|'
                 r'hydroxychloroquine blood|vitamin|exercise|diet|education|depress', re.I)
root = ET.parse(os.path.join(D, 'raw', 'pubmed_efetch.xml')).getroot()
out = []
for a in root.findall('.//PubmedArticle'):
    pmid = a.findtext('.//MedlineCitation/PMID')
    title = ''.join(a.find('.//ArticleTitle').itertext())
    yr = a.findtext('.//JournalIssue/PubDate/Year') or (a.findtext('.//JournalIssue/PubDate/MedlineDate') or '')[:4]
    ncts = sorted({x.text for x in a.findall('.//DataBank/AccessionNumberList/AccessionNumber') if x.text and x.text.startswith('NCT')})
    ncts = sorted(set(ncts) | set(re.findall(r'NCT\d{8}', title)))
    status = 'in60' if any(n in known for n in ncts) else ('secondary_title' if SEC.search(title) else 'screen')
    out.append(dict(pmid=pmid, year=yr, ncts=';'.join(ncts), status=status, title=title))
with open(os.path.join(D, 'pubmed_screen.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
from collections import Counter; print(Counter(r['status'] for r in out))
