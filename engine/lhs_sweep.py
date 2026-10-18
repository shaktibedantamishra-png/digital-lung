# Joint sensitivity (not pre-specified; PD7): 500 Latin-hypercube draws over the OAT ranges.
import sys, numpy as np, pandas as pd, time
sys.path.insert(0, 'engine'); import digital_lung as dl
from scipy.stats import qmc
R = dict(r=(0.1, 0.5), sigma=(3, 8), kappa=(0.33, 0.95), R_aw=(5, 15), VD_kg=(1.5, 3.0), eps_star=(1.5, 2.0))
X = qmc.LatinHypercube(d=len(R), seed=20261007).random(500)
lo = np.array([v[0] for v in R.values()]); hi = np.array([v[1] for v in R.values()])
X = lo + X * (hi - lo)
pop = dl.resample(dl.load_pairs()); G = dl.grid(); out = []; t0 = time.time()
for i, x in enumerate(X):
    p = dict(dl.BASE); p.update(dict(zip(R, x)))
    s = dl.pop_summary(dl.evaluate(pop, G, p))
    mins = [dl.min_cost(s, wd, wo) for wd in dl.W_GRID for wo in dl.W_GRID]
    floor = sum((m.vt_kg == 4) and (m.peep == 5) and (m.rr == 15) for m in mins)
    Mg = dl.matched(kind='va', p=p); Mg = Mg[Mg.matched]
    sv = dl.pop_summary(dl.evaluate(pop, Mg[['vt_kg', 'peep', 'rr']], p))
    m0 = sv.loc[sv.mp.idxmin()]; m11 = dl.min_cost(sv, 1, 1)
    out.append(dict(draw=i, **dict(zip(R, x)), n_w_floor=floor, va_mp_vt=m0.vt_kg, va_mp_peep=m0.peep, va_w11_vt=m11.vt_kg, va_w11_peep=m11.peep))
    if i % 50 == 0: print(i, round(time.time() - t0), flush=True)
pd.DataFrame(out).to_csv('lhs_joint_sensitivity.csv', index=False); print('done')
