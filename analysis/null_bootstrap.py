"""Noise-only null, PEEP-error variants, stay-bootstrap CIs for the excess over the null (PD4/PD6),
and set-PEEP hysteresis incl. the <=5 min sensitivity (PD8). Seed 20261007. Run from the repository root."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import hysteresis as hy
SEED = 20261007
OUT = os.path.join(os.environ.get('REPO_ROOT', os.getcwd()), 'results_reproduced'); os.makedirs(OUT, exist_ok=True)

def prep(x, phen):
    x = ph.vt_bin(x, 'pbw'); x = x[x.stay_id.isin(phen.stay_id)].copy(); x['pl'] = x.peep.round()
    a = x.groupby(['stay_id', 'vbin']).agg(nlev=('pl', 'nunique'), lo=('pl', 'min'), hi=('pl', 'max'))
    q = a[(a.nlev >= 2) & ((a.hi - a.lo) >= 3)].reset_index(); q['gid'] = np.arange(len(q))
    C = x.groupby('stay_id').crs.median()
    x = x.merge(q[['stay_id', 'vbin', 'gid']], on=['stay_id', 'vbin']); x['C'] = x.stay_id.map(C)
    return x, q

def slopes(g, X, Y, G):
    n = np.bincount(g, minlength=G); sx = np.bincount(g, X, G); sy = np.bincount(g, Y, G)
    sxx = np.bincount(g, X*X, G); sxy = np.bincount(g, X*Y, G)
    den = n*sxx - sx**2
    return (n*sxy - sx*sy) / np.where(den > 0, den, np.nan)

def psd(x):
    g = x.groupby(['stay_id', 'vbin', 'pl']).pplat; v = g.var(ddof=1); n = g.size(); ok = n > 1
    return float(np.sqrt(((n[ok]-1)*v[ok]).sum() / (n[ok]-1).sum()))

def excess_full(x, phen, db, B=200, NB=2000, seed=SEED):
    x, q = prep(x, phen); G = len(q); g = x.gid.values; pe = x.peep.values
    base = pe + x.vt.values / x.C.values; se = psd(x); sid = q.stay_id.values
    cls = lambda s: pd.Series(s).groupby(sid).median()
    obs = cls(slopes(g, pe, x.pplat.values, G))
    rg = np.random.default_rng(seed); pn = np.zeros(len(obs)); po = np.zeros(len(obs)); R = []; M = []
    for b in range(B):
        m = cls(slopes(g, pe, base + rg.normal(0, se, len(x)), G)).values
        pn += (m < 0.9); po += (m > 1.1); R.append((m < 0.9).mean()); M.append(np.median(m))
    pn /= B; po /= B
    o = (obs.values < 0.9).astype(float); d = o - pn
    bs = np.array([d[rg.integers(0, len(d), len(d))].mean() for _ in range(NB)]) * 100
    k = int(o.sum()); lo, hi = ph.wilson(k, len(o)); pc = lambda a, p: 100*np.percentile(a, p)
    return dict(db=db, n=len(obs), s_e=round(se, 2), obs_rec=round(100*o.mean(), 1),
                obs_rec_ci=f"{100*lo:.1f} to {100*hi:.1f}", null_rec=round(100*pn.mean(), 1),
                null_rec_pi=f"{pc(R, 2.5):.1f} to {pc(R, 97.5):.1f}", excess=round(100*d.mean(), 1),
                excess_ci=f"{np.percentile(bs, 2.5):.1f} to {np.percentile(bs, 97.5):.1f}",
                obs_over=round(100*(obs.values > 1.1).mean(), 1), null_over=round(100*po.mean(), 1),
                obs_median=round(float(obs.median()), 2), null_median=round(float(np.mean(M)), 2))

if __name__ == '__main__':
    ev, prim, per, SET, NOI = load_sets()
    rows = []
    for db in ['mimic', 'hirid']:
        rows.append(dict(analysis='Set PEEP, all eligible stays', **excess_full(SET[db], ph.phenotype(SET[db], db), db)))
        rows.append(dict(analysis='Set PEEP, no-intrinsic-PEEP subset', **excess_full(NOI[db], ph.phenotype(NOI[db], db), db)))
    rows.append(dict(analysis='Set PEEP (eICU-CRD primary)', **excess_full(prim['eicu'], per['eicu'], 'eicu')))
    P6 = pd.DataFrame(rows); P6.to_csv(os.path.join(OUT, 'pd6_excess_full.csv'), index=False)
    P6[['db', 'analysis', 'n', 's_e', 'obs_rec', 'null_rec', 'excess', 'excess_ci']].rename(
        columns={'obs_rec': 'obs', 'null_rec': 'null', 'excess_ci': 'ci'}).to_csv(os.path.join(OUT, 'excess_bootstrap_ci.csv'), index=False)
    H = []
    cases = [('Measured total PEEP (primary)', lambda db: prim[db]), ('Set PEEP', lambda db: SET[db]),
             ('Set PEEP, set-PEEP lag <=5 min (PD8)', lambda db: SET[db][SET[db].lag_peep_set <= 5]),
             ('Set PEEP, no-intrinsic-PEEP subset (exploratory)', lambda db: NOI[db])]
    for db in ['mimic', 'hirid']:
        for lab, f in cases:
            o, _ = hy.paired(hy.steps(f(db)[COLS]), db)
            lo, hi = [float(v) for v in o['ci'].split(' to ')]
            H.append(dict(db=db, analysis=lab, n_paired=o['n_paired'], median_increase=o['median_increase'],
                          median_decrease=o['median_decrease'], diff_increase_minus_decrease=round(-o['median_diff'], 2),
                          ci_increase_minus_decrease=f"{-hi:.2f} to {-lo:.2f}", wilcoxon_p=o['wilcoxon_p']))
    pd.DataFrame(H).to_csv(os.path.join(OUT, 'hysteresis_setpeep_pd6_pd8.csv'), index=False)
    print('done', OUT)
