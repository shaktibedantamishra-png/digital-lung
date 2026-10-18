"""Direction asymmetry of PEEP steps (analysis_spec_v1 S7)."""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from phenotype import vt_bin, SEED

def steps(d, binning="pbw", min_step=2.0, max_gap_min=240.0):
    d = vt_bin(d, binning).sort_values(["stay_id", "t"])
    g = d.groupby("stay_id")
    s = pd.DataFrame(dict(stay_id=d.stay_id, dpeep=g.peep.diff(), dpplat=g.pplat.diff(),
                          gap=g.t.diff().dt.total_seconds() / 60 if np.issubdtype(d.t.dtype, np.datetime64) else g.t.diff(),
                          same_bin=g.vbin.diff() == 0))
    s = s[s.same_bin & (s.dpeep.abs() >= min_step) & (s.gap <= max_gap_min)].copy()
    s["slope"] = s.dpplat / s.dpeep
    s["direction"] = np.where(s.dpeep > 0, "increase", "decrease")
    return s

def paired(s, db, B=2000, seed=SEED):
    m = s.groupby(["stay_id", "direction"]).slope.median().unstack()
    n_inc, n_dec = int(m.get("increase").notna().sum()), int(m.get("decrease").notna().sum())
    p = m.dropna(subset=["increase", "decrease"])
    diff = (p.decrease - p.increase).values
    out = dict(db=db, n_steps=len(s), stays_increase=n_inc, stays_decrease=n_dec, n_paired=len(diff))
    if len(diff) < 10: return out, p
    rng = np.random.default_rng(seed)
    bs = np.array([np.median(diff[rng.integers(0, len(diff), len(diff))]) for _ in range(B)])
    out.update(median_increase=round(float(p.increase.median()), 2), median_decrease=round(float(p.decrease.median()), 2),
               median_diff=round(float(np.median(diff)), 2), ci=f"{np.percentile(bs,2.5):.2f} to {np.percentile(bs,97.5):.2f}",
               wilcoxon_p=f"{wilcoxon(diff).pvalue:.1e}", pct_decrease_lower=round(100 * float((diff < 0).mean()), 1))
    return out, p
