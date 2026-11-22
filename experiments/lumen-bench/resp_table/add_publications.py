"""Append round-90 published week-104 rows (not in the registry) -> resp_table_plus_pub.csv.
NOBILITY: Furie 2022 Ann Rheum Dis, 26/63 vs 14/62 (counts). TULIP-LN: Jayne 2023 Lupus Sci Med, 27.3/18.6/17.8%
-> counts rounded on randomized N (51/45/49); approximate, year-2 population selected (see round 90)."""
import pandas as pd
R=pd.read_csv('resp_table.csv')
P=pd.DataFrame([
 dict(trial='NOBILITY',nct='NCT02550652',arm='obinutuzumab',molecule='obinutuzumab',target='CD20',cls='slow',week=104,r_drug=26,n_drug=63,r_ctrl=14,n_ctrl=62,source='publication:Furie2022'),
 dict(trial='TULIP-LN1',nct='NCT02547922',arm='anifro_int',molecule='anifrolumab',target='IFNAR',cls='slow',week=104,r_drug=14,n_drug=51,r_ctrl=9,n_ctrl=49,source='publication:Jayne2023(pct)'),
 dict(trial='TULIP-LN1',nct='NCT02547922',arm='anifro_basic',molecule='anifrolumab',target='IFNAR',cls='slow',week=104,r_drug=8,n_drug=45,r_ctrl=9,n_ctrl=49,source='publication:Jayne2023(pct)')])
pd.concat([R,P],ignore_index=True).to_csv('resp_table_plus_pub.csv',index=False);print(len(R)+len(P),'rows')
