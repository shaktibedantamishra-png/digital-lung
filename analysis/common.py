"""Shared setup: event tables, primary and set-PEEP analysis sets. Run from the repository root.
Requires pipe_out/events_{db}.parquet, produced by pipeline/build_events.py."""
import os, sys
import pandas as pd
ROOT = os.environ.get('REPO_ROOT', os.getcwd())
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import filters, phenotype as ph, build_events as be
DBO = ['mimic', 'eicu', 'hirid']
LAB = {'mimic': 'MIMIC-IV', 'eicu': 'eICU-CRD', 'hirid': 'HiRID'}
COLS = ['stay_id', 't', 'pplat', 'peep', 'vt', 'pbw', 'mode', 'crs']

def apply_ext(e, db, tol=60.0, vt='vt_obs'):
    d = e[(e.lag_peep <= tol)]; d = d[d[f'lag_{vt}'] <= tol]
    d = d.assign(vt=d[vt].astype('float64'), peep=d.peep.astype('float64'), pplat=d.plateau_pressure.astype('float64'))
    d = d[d.pplat.between(5, 50) & d.peep.between(0, 30) & d.vt.between(150, 1200) & (d.pplat > d.peep)]
    d = d.assign(mode=d.mode_class.where(d.lag_mode_class <= filters.MODE_TOL, 'unknown').fillna('unknown'))
    d = d[~d['mode'].isin(['spontaneous', 'not_invasive'])]
    return d.assign(crs=d.vt / (d.pplat - d.peep)).reset_index(drop=True)

def load_sets():
    ev = {db: pd.read_parquet(os.path.join(ROOT, 'pipe_out', f'events_{db}.parquet')) for db in DBO}
    prim = {db: filters.apply_spec(ev[db], db, log=[]) for db in DBO}
    per = {db: ph.phenotype(prim[db], db) for db in DBO}
    SET, NOI = {}, {}
    for db in ['mimic', 'hirid']:
        x = apply_ext(ev[db], db)
        y = be.asof(x.copy(), be.load(db, 'peep_set', set(x.stay_id)), 'peep_set')
        y['has'] = (y.lag_peep_set <= 60) & y.peep_set.between(0, 30)
        y['agree'] = y.has & ((y.peep - y.peep_set).abs() <= 0.5)
        st = y.groupby('stay_id').agg(a=('has', 'all'), b=('agree', 'all'))
        keep = set(st.index[st.a & st.b])
        ys = y[y.has].assign(peep_tot=lambda d: d.peep, peep=lambda d: d.peep_set.astype(float))
        ys = ys[ys.pplat > ys.peep].copy(); ys['crs'] = ys.vt / (ys.pplat - ys.peep)
        SET[db] = ys.reset_index(drop=True)
        NOI[db] = ys[ys.stay_id.isin(keep)].reset_index(drop=True)
    return ev, prim, per, SET, NOI
