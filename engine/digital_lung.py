"""Digital-lung trial engine, rebuilt to engine_spec_v1.md (locked, option A).

Units: volumes in L, pressures in cmH2O, compliance in L/cmH2O, RR in /min, power in J/min.
All equations in this file are the equations reported in Additional file 1, section S1.
"""
import numpy as np, pandas as pd
from scipy.stats import norm

BASE = dict(EL_spec=13.5, kappa=0.7, r=0.3, mu_c=5.0, g=15.0, sigma=5.0,
            eps_star=1.5, R_aw=10.0, IE=0.5, VD_kg=2.2)
SEED = 20260929
VT_GRID = [4, 6, 8, 10]; PEEP_GRID = [5, 8, 10, 12, 15, 18, 20]; RR_GRID = [15, 20, 25]
W_GRID = np.round(np.arange(0, 4.01, 0.5), 2)
DBS = ['mimic', 'eicu', 'hirid']


# ---------- E1 population ----------
def load_pairs(crs_file='compliance_perstay.parquet', pbw_dir='pipe_out'):
    c = pd.read_parquet(crs_file)
    out = []
    for db in DBS:
        p = pd.read_parquet(f'{pbw_dir}/cohort_{db}_pbw.parquet')[['stay_id', 'pbw']]
        x = c[c.db == db].merge(p, on='stay_id', how='inner').dropna(subset=['crs', 'pbw'])
        out.append(x)
    return pd.concat(out, ignore_index=True)


def resample(pairs, n=2000, seed=SEED, method='empirical'):
    rng = np.random.default_rng(seed); out = []
    for i, db in enumerate(DBS):
        x = pairs[pairs.db == db]
        if method == 'empirical':
            s = x.sample(n=n, replace=True, random_state=rng.integers(2**31)).reset_index(drop=True)
            s = s[['db', 'crs', 'pbw']]
        elif method == 'copula':
            u = x[['crs', 'pbw']].rank(pct=True).clip(1e-6, 1 - 1e-6)
            z = norm.ppf(u.values); rho = np.corrcoef(z.T)[0, 1]
            zz = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=n)
            uu = norm.cdf(zz)
            s = pd.DataFrame({'db': db, 'crs': np.quantile(x.crs, uu[:, 0]),
                              'pbw': np.quantile(x.pbw, uu[:, 1])})
        out.append(s)
    pop = pd.concat(out, ignore_index=True)
    pop['C'] = pop.crs / 1000.0          # mL/cmH2O -> L/cmH2O
    return pop


# ---------- E2 lung state ----------
def frc(C, p):  return p['EL_spec'] * C / p['kappa']


def fractions(peep, pplat, p):
    hi_open = np.minimum(peep, pplat - p['g'])
    f_open = norm.cdf((hi_open - p['mu_c']) / p['sigma'])
    f_cyc = np.clip(norm.cdf((pplat - p['g'] - p['mu_c']) / p['sigma'])
                    - norm.cdf((peep - p['mu_c']) / p['sigma']), 0, None)
    return f_open, f_cyc, 1 - f_open - f_cyc


# ---------- E3 mechanical power ----------
def mp(C, vt, rr, peep, p):
    K = (1 + p['IE']) / (60 * p['IE'])
    return 0.098 * rr * (0.5 * vt**2 / C + rr * K * p['R_aw'] * vt**2 + peep * vt)


# ---------- E4 penalties ----------
def penalties(C, vt, rr, peep, p):
    pplat = peep + vt / C
    f_open, f_cyc, _ = fractions(peep, pplat, p)
    F = frc(C, p); VR = p['r'] * F
    p_dr = 0.098 * rr * VR * f_cyc * (p['mu_c'] + p['g'])
    A = F + VR * f_open; Vstar = p['eps_star'] * A
    Vpeep = C * peep; Vei = Vpeep + vt; EL = p['kappa'] / C
    p_od = np.where(Vei > Vstar, 0.098 * rr * 0.5 * EL * (Vei**2 - np.maximum(Vpeep, Vstar)**2), 0.0)
    strain = Vei / A
    return p_dr, p_od, dict(pplat=pplat, f_open=f_open, f_cyc=f_cyc, strain=strain)


def evaluate(pop, strategies, p=BASE):
    """strategies: DataFrame with vt_kg, peep, rr. Returns long per-patient x strategy frame."""
    C = pop.C.values[:, None]; pbw = pop.pbw.values[:, None]
    vt = strategies.vt_kg.values[None, :] * pbw / 1000.0
    peep = strategies.peep.values[None, :].astype(float); rr = strategies.rr.values[None, :].astype(float)
    m = mp(C, vt, rr, peep, p); dr, od, aux = penalties(C, vt, rr, peep, p)
    n, k = m.shape
    df = pd.DataFrame({'pid': np.repeat(np.arange(n), k), 'db': np.repeat(pop.db.values, k),
                       'sid': np.tile(np.arange(k), n)})
    for c in ['vt_kg', 'peep', 'rr']: df[c] = np.tile(strategies[c].values, n)
    df['mp'] = m.ravel(); df['p_dr'] = dr.ravel(); df['p_od'] = od.ravel()
    df['pplat'] = aux['pplat'].ravel(); df['dp'] = (vt / C).ravel(); df['strain'] = aux['strain'].ravel()
    df['f_cyc'] = aux['f_cyc'].ravel()
    return df


def grid(vts=VT_GRID, peeps=PEEP_GRID, rrs=RR_GRID):
    return pd.DataFrame([(v, pe, r) for v in vts for pe in peeps for r in rrs], columns=['vt_kg', 'peep', 'rr'])


def matched(vts=VT_GRID, peeps=PEEP_GRID, kind='mv', p=BASE, VE_kg=100.0, rr_lo=10, rr_hi=45):
    rows = []
    for v in vts:
        if kind == 'mv': rr = VE_kg / v
        else:
            VA_kg = 20 * (6 - p['VD_kg']); rr = VA_kg / (v - p['VD_kg']) if v > p['VD_kg'] else np.inf
        rc = float(np.clip(rr, rr_lo, rr_hi))
        for pe in peeps: rows.append((v, pe, rc, bool(np.isclose(rc, rr))))
    return pd.DataFrame(rows, columns=['vt_kg', 'peep', 'rr', 'matched'])


def pop_summary(ev, by=('vt_kg', 'peep', 'rr')):
    return ev.groupby(list(by), as_index=False)[['mp', 'p_dr', 'p_od', 'pplat', 'dp', 'strain']].mean()


def min_cost(summ, w_dr, w_od):
    J = summ.mp + w_dr * summ.p_dr + w_od * summ.p_od
    r = summ.loc[J.idxmin()].copy(); r['J'] = J.min(); return r
