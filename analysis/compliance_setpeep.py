"""Static compliance recomputed with set PEEP (mimic, hirid) against total PEEP on the same stays (PD6).
Run from the repository root."""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
OUT = os.path.join(os.environ.get('REPO_ROOT', os.getcwd()), 'results_reproduced'); os.makedirs(OUT, exist_ok=True)
if __name__ == '__main__':
    ev, prim, per, SET, NOI = load_sets()
    rows = []
    for db in ['mimic', 'hirid']:
        d = SET[db].assign(crs_tot=lambda d: d.vt / (d.pplat - d.peep_tot))
        s = d.groupby('stay_id').agg(set_=('crs', 'median'), tot=('crs_tot', 'median'))
        dd = s.set_ - s.tot
        rows.append(dict(db=db, stays=len(s), median_crs_total_all=round(prim[db].groupby('stay_id').crs.median().median(), 1),
                         median_crs_total=round(s.tot.median(), 1), median_crs_set=round(s.set_.median(), 1),
                         median_diff=round(dd.median(), 2), iqr=f"{dd.quantile(.25):.2f} to {dd.quantile(.75):.2f}",
                         median_rel_pct=round((100*(s.set_/s.tot - 1)).median(), 1),
                         pct_stays_set_lower=round(100*(dd < -0.5).mean(), 1),
                         mean_tot_minus_set_peep=round((d.peep_tot - d.peep).mean(), 2)))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'compliance_setpeep.csv'), index=False)
    print('done', OUT)
