"""Rebuild the primary-endpoint arm table (log OR, SE) used by lumen_forecast.py.
resp_table.csv was never committed (*.csv is git-ignored) and is absent from this workspace, so the
primary rows are re-entered from the logged round-86 output (result86.txt, 2 dp) and the round-90
scored week-104 rows (r90_scored.txt). Counts are not needed for a prior-only forecast."""
import pandas as pd
T={'abatacept':'CD80/86','iscalimab':'CD40','zetomipzomib':'immunoproteasome','ocrelizumab':'CD20',
 'voclosporin':'calcineurin','ALXN2050':'factor D','ravulizumab':'C5','guselkumab':'IL-23','anifrolumab':'IFNAR',
 'secukinumab':'IL-17A','obinutuzumab':'CD20','rituximab':'CD20','BIIB023':'TWEAK','BI655064':'CD40','belimumab':'BAFF'}
rows="""ACCESS abatacept abatacept 24 0.11 0.37
CFZ533 iscalimab iscalimab 24 0.57 0.64
PALIZADE zeto_60 zetomipzomib 25 0.33 0.82
PALIZADE zeto_30 zetomipzomib 25 -1.06 1.05
BELONG ocre_1000 ocrelizumab 48 -0.14 0.35
BELONG ocre_400 ocrelizumab 48 0.33 0.33
AURA-LV vocl_low voclosporin 48 1.12 0.33
AURA-LV vocl_high voclosporin 48 0.73 0.33
ALXN2050 fD_120 ALXN2050 50 -2.03 1.60
ALXN2050 fD_180 ALXN2050 50 -1.02 0.84
SANCTUARY ravulizumab ravulizumab 50 0.42 0.58
ORCHID-LN guselkumab guselkumab 52 -0.95 1.07
TULIP-LN1 anifro_int anifrolumab 52 0.61 0.41
SELUNE secukinumab secukinumab 52 -0.55 0.32
BMS-LN abata_30 abatacept 52 0.13 0.49
NOBILITY obinutuzumab obinutuzumab 52 0.60 0.40
AURORA voclosporin voclosporin 52 0.86 0.23
BMS-LN abata_10 abatacept 52 0.35 0.48
LUNAR rituximab rituximab 52 -0.20 0.37
ATLAS biib_20 BIIB023 52 0.27 0.74
ATLAS biib_3 BIIB023 52 0.25 0.74
BI655064 bi_120 BI655064 52 -0.36 0.54
TULIP-LN1 anifro_basic anifrolumab 52 -0.84 0.50
BI655064 bi_240 BI655064 52 -0.10 0.44
BI655064 bi_180 BI655064 52 -0.09 0.54
REGENCY obinutuzumab obinutuzumab 76 0.57 0.25
CALIBRATE belimumab belimumab 96 0.45 0.62
ALLURE abatacept abatacept 104 0.29 0.20
BLISS-LN belimumab belimumab 104 0.55 0.22
NOBILITY obinutuzumab obinutuzumab 104 0.86 0.39
TULIP-LN1 anifro_int anifrolumab 104 0.50 0.47
TULIP-LN1 anifro_basic anifrolumab 104 -0.03 0.52"""
R=pd.DataFrame([r.split() for r in rows.splitlines()],columns=['trial','arm','molecule','week','lor','se'])
R[['week','lor','se']]=R[['week','lor','se']].astype(float);R['target']=R.molecule.map(T)
R.to_csv('arm_table_reconstructed.csv',index=False);print(len(R),'rows')
