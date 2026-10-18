"""Apply analysis_spec_v1 S2-S3 to plateau-anchored events, with an attrition log."""
import numpy as np, pandas as pd

MODE_TOL = 720.0  # D3

def apply_spec(e, db, tol=60.0, vt="vt_obs", passivity="permissive", log=None):
    att = []
    def step(name, d): att.append(dict(db=db, step=name, rows=len(d), stays=d.stay_id.nunique())); return d
    d = step("plateau rows (cohort)", e)
    d = step(f"PEEP within {tol:g} min", d[d.lag_peep <= tol])
    d = step(f"{vt} within {tol:g} min", d[d[f"lag_{vt}"] <= tol])
    d = d.assign(vt=d[vt].astype("float64"), peep=d.peep.astype("float64"), pplat=d.plateau_pressure.astype("float64"))
    ok = d.pplat.between(5, 50) & d.peep.between(0, 30) & d.vt.between(150, 1200) & (d.pplat > d.peep)
    d = step("plausibility bounds", d[ok])
    mode = d.mode_class.where(d.lag_mode_class <= MODE_TOL, "unknown").fillna("unknown")
    d = d.assign(mode=mode)
    if passivity == "permissive":
        d = step("exclude spontaneous / non-invasive mode", d[~d["mode"].isin(["spontaneous", "not_invasive"])])
    elif passivity == "strict":
        d = step("controlled mode only", d[d["mode"] == "controlled"])
    elif passivity == "rr":
        d = d[~d["mode"].isin(["spontaneous", "not_invasive"])]
        rr_ok = (d.lag_rr_tot <= tol) & (d.lag_rr_set <= tol) & (d.rr_tot <= d.rr_set + 2)
        d = step("permissive + RR total <= RR set + 2", d[rr_ok])
    d = d.assign(crs=d.vt / (d.pplat - d.peep))
    if log is not None: log.extend(att)
    return d[["stay_id", "t", "pplat", "peep", "vt", "pbw", "mode", "crs"]].reset_index(drop=True)
