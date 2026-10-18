"""PEEP-response phenotype (analysis_spec_v1 S6, D6) and shared statistics (S8)."""
import numpy as np, pandas as pd
from scipy.stats import chi2_contingency

SEED = 20260925
CLASSES = ["recruitment-consistent", "neutral", "overdistension-consistent"]

def vt_bin(d, binning="pbw"):
    if binning == "pbw":
        d = d[d.pbw.notna()]
        return d.assign(vbin=np.floor(d.vt / d.pbw).astype(int))
    return d.assign(vbin=np.floor(d.vt / 50.0).astype(int))

def classify(s):
    return np.select([s < 0.9, s <= 1.1], CLASSES[:2], CLASSES[2])

def phenotype(d, db, binning="pbw", min_span=3.0):
    d = vt_bin(d, binning).assign(peep_lvl=lambda x: x.peep.round())
    t0 = d.groupby("stay_id").t.min().rename("t0")          # first analysable plateau record
    g = d.groupby(["stay_id", "vbin"])
    agg = g.agg(n=("pplat", "size"), nlev=("peep_lvl", "nunique"), lo=("peep_lvl", "min"), hi=("peep_lvl", "max"),
                sx=("peep", "sum"), sy=("pplat", "sum"), t_first=("t", "min"))
    d = d.assign(xx=d.peep**2, xy=d.peep * d.pplat)
    s2 = d.groupby(["stay_id", "vbin"])[["xx", "xy"]].sum()
    agg = agg.join(s2)
    q = agg[(agg.nlev >= 2) & ((agg.hi - agg.lo) >= min_span)].copy()
    den = q.n * q.xx - q.sx**2
    q["slope"] = (q.n * q.xy - q.sx * q.sy) / den
    q = q[np.isfinite(q.slope) & (den > 0)].reset_index()
    per = q.groupby("stay_id").agg(slope=("slope", "median"), n_bins=("slope", "size"), t_first=("t_first", "min")).reset_index()
    per = per.merge(t0, on="stay_id")
    h = (per.t_first - per.t0) / 60.0
    per["time_window"] = pd.cut(h, [-np.inf, 2, 6, 24, np.inf], right=False, labels=["<2 h", "2-6 h", "6-24 h", ">24 h"])
    per["phenotype"] = classify(per.slope.values)
    per["db"] = db
    return per[["stay_id", "db", "slope", "phenotype", "n_bins", "time_window"]]

def wilson(k, n, z=1.959964):
    p = k / n; den = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / den; h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return c - h, c + h

def boot_median_ci(x, B=2000, seed=SEED):
    rng = np.random.default_rng(seed); x = np.asarray(x); n = len(x)
    m = np.empty(B)
    for b in range(B): m[b] = np.median(x[rng.integers(0, n, n)])
    return np.percentile(m, [2.5, 97.5])

def cramers_v(tab):
    chi2, p, dof, _ = chi2_contingency(tab)
    n = tab.values.sum(); k = min(tab.shape) - 1
    return chi2, p, dof, np.sqrt(chi2 / (n * k))

def cramers_v_boot(per, B=2000, seed=SEED):
    rng = np.random.default_rng(seed); vs = []
    groups = {db: g.phenotype.values for db, g in per.groupby("db")}
    for b in range(B):
        rows = [pd.Series(v[rng.integers(0, len(v), len(v))]).value_counts().reindex(CLASSES, fill_value=0).rename(db)
                for db, v in groups.items()]
        vs.append(cramers_v(pd.concat(rows, axis=1))[3])
    return np.percentile(vs, [2.5, 97.5])

def table3(per):
    rows = []
    for db, g in per.groupby("db", sort=False):
        n = len(g); r = dict(db=db, n=n)
        for c in CLASSES:
            k = int((g.phenotype == c).sum()); lo, hi = wilson(k, n)
            r[f"{c}_k"] = k; r[f"{c}_pct"] = round(100 * k / n, 1); r[f"{c}_ci"] = f"{100*lo:.1f} to {100*hi:.1f}"
        mlo, mhi = boot_median_ci(g.slope.values)
        r["median_slope"] = round(float(np.median(g.slope)), 2); r["median_slope_ci"] = f"{mlo:.2f} to {mhi:.2f}"
        r["iqr_slope"] = f"{g.slope.quantile(.25):.2f} to {g.slope.quantile(.75):.2f}"
        rows.append(r)
    return pd.DataFrame(rows)
