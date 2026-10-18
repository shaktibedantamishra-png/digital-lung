"""Runs every analysis in engine_spec_v1.md E5 (+E9). Outputs to trial_out/."""
import sys, os, json, hashlib, itertools, numpy as np, pandas as pd
sys.path.insert(0, 'engine'); import digital_lung as dl
O = 'trial_out'; os.makedirs(O, exist_ok=True)
sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
log = {'spec_sha256': sha('engine_spec_v1.md'), 'engine_sha256': sha('engine/digital_lung.py'),
       'tests_sha256': sha('tests/test_engine.py'), 'seed': dl.SEED}
W = list(itertools.product(dl.W_GRID, dl.W_GRID))

pairs = dl.load_pairs(); log['pairs_per_db'] = pairs.groupby('db').size().to_dict()
pop = dl.resample(pairs)
log['pop_crs_median'] = pop.groupby('db').crs.median().round(2).to_dict()
log['pop_pbw_median'] = pop.groupby('db').pbw.median().round(2).to_dict()
G = dl.grid()


def summaries(ev):
    s = dl.pop_summary(ev).assign(db='pooled')
    return pd.concat([s] + [dl.pop_summary(ev[ev.db == d]).assign(db=d) for d in dl.DBS], ignore_index=True)


def sweep(summ, label, by_db=True):
    rows = []
    for d, s in summ.groupby('db'):
        if not by_db and d != 'pooled': continue
        s = s.reset_index(drop=True)
        for wd, wo in W:
            m = dl.min_cost(s, wd, wo)
            rows.append(dict(analysis=label, db=d, w_dr=wd, w_od=wo, vt_kg=m.vt_kg, peep=m.peep, rr=m.rr, J=m.J,
                             mp=m.mp, p_dr=m.p_dr, p_od=m.p_od))
    return pd.DataFrame(rows)


# E5.1 main grid
ev = dl.evaluate(pop, G); S = summaries(ev); S.to_csv(f'{O}/trial_main.csv', index=False)
sw = sweep(S, 'main'); sw.to_csv(f'{O}/weight_sweep.csv', index=False)

# E5.2 / E5.3 matched ventilation
mres = []
for kind in ['mv', 'va']:
    Mg = dl.matched(kind=kind); e = dl.evaluate(pop, Mg[['vt_kg', 'peep', 'rr']])
    e = e.merge(Mg.reset_index().rename(columns={'index': 'sid'})[['sid', 'matched']], on='sid')
    s = summaries(e).merge(Mg, on=['vt_kg', 'peep', 'rr'])
    s.to_csv(f'{O}/trial_matched_{kind}.csv', index=False)
    for scope, ss in [('peep10', s[s.peep == 10]), ('allpeep', s)]:
        for only in [True, False]:
            x = ss[ss.matched] if only else ss
            r = sweep(x, f'matched_{kind}_{scope}_{"matchedonly" if only else "all"}', by_db=False); mres.append(r)
pd.concat(mres).to_csv(f'{O}/matched_sweep.csv', index=False)

# E5.4 protective limits, per patient
feas = ev[(ev.pplat <= 30) & (ev.dp <= 15)].copy()
log['protective_patients_with_any_feasible'] = int(feas.pid.nunique()); log['n_patients'] = int(len(pop))
prot = []
for wd, wo in [(0, 0), (1, 0), (0, 1), (1, 1)]:
    J = feas.mp + wd * feas.p_dr + wo * feas.p_od
    best = feas.loc[J.groupby(feas.pid).idxmin()]
    for c in ['vt_kg', 'peep', 'rr']:
        t = best.groupby(['db', c]).size().rename('n').reset_index().rename(columns={c: 'value'})
        t['setting'] = c; t['w_dr'] = wd; t['w_od'] = wo; prot.append(t)
pd.concat(prot).to_csv(f'{O}/protective_perpatient.csv', index=False)

# E5.5 one-at-a-time sensitivity
variants = [('base', {})] + [(f'{k}={v}', {k: v}) for k, vs in
            dict(r=[0.1, 0.5], sigma=[3, 8], kappa=[0.33, 0.95], eps_star=[2.0], R_aw=[5, 15], VD_kg=[1.5, 3.0]).items()
            for v in vs] + [('copula', 'copula')]
oat = []
for name, ch in variants:
    p = dict(dl.BASE); pp = pop
    if ch == 'copula': pp = dl.resample(pairs, method='copula')
    else: p.update(ch)
    s = dl.pop_summary(dl.evaluate(pp, G, p)).assign(db='pooled')
    m11 = dl.min_cost(s, 1, 1); swp = sweep(s, name, by_db=False)
    row = dict(variant=name, vt_kg_w11=m11.vt_kg, peep_w11=m11.peep, rr_w11=m11.rr,
               n_w_peep_min=int((swp.peep == 5).sum()), n_w_peep_max=int((swp.peep == 20).sum()),
               n_w_vt_min=int((swp.vt_kg == 4).sum()), peep_w00=dl.min_cost(s, 0, 0).peep)
    for kind in ['mv', 'va']:
        Mg = dl.matched(kind=kind, p=p); Mg = Mg[Mg.matched]
        sm = dl.pop_summary(dl.evaluate(pp, Mg[['vt_kg', 'peep', 'rr']], p))
        row[f'{kind}_vt_w11_allpeep'] = dl.min_cost(sm, 1, 1).vt_kg
        row[f'{kind}_vt_w00_allpeep'] = dl.min_cost(sm, 0, 0).vt_kg
    oat.append(row)
pd.DataFrame(oat).to_csv(f'{O}/sensitivity_oat.csv', index=False)

# bootstrap of the population
bt = []
for b in range(200):
    pb = dl.resample(pairs, seed=dl.SEED + 1 + b)
    s = dl.pop_summary(dl.evaluate(pb, G)).assign(db='pooled')
    for wd, wo in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        m = dl.min_cost(s, wd, wo); bt.append(dict(b=b, w_dr=wd, w_od=wo, vt_kg=m.vt_kg, peep=m.peep, rr=m.rr))
bt = pd.DataFrame(bt)
bt.groupby(['w_dr', 'w_od', 'vt_kg', 'peep', 'rr']).size().rename('n_of_200').reset_index().to_csv(
    f'{O}/bootstrap_freq.csv', index=False)
json.dump(log, open(f'{O}/run_log.json', 'w'), indent=1, default=str)
print(json.dumps(log, default=str))
