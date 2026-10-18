"""Plateau-anchored event construction (analysis_spec_v1.md S1-S4).
Run: python pipeline/build_events.py  -> events_<db>.parquet, cohort_<db>_pbw.parquet, attrition_events.csv
Events keep lag (min) for every joined variable; tolerance is applied downstream so caps can be varied.
"""
import os, numpy as np, pandas as pd, pyarrow.parquet as pq, pyarrow.dataset as ds
P = os.environ.get("SUPERICU_PARQUET", os.environ['SUPERICU_DATA'])
OUT = os.environ.get("PIPE_OUT", "pipe_out"); os.makedirs(OUT, exist_ok=True)
DBS = ["mimic", "eicu", "hirid"]
MODE_MAP = pd.read_csv(os.environ.get("MODE_MAP", "mode_mapping.csv"))

def tmin(df, db):
    """Return time in float minutes (eICU offsets already minutes)."""
    if "charttime" in df: return df["charttime"].astype("datetime64[us]").astype("int64") / 6e7
    if "charttime_offset" in df: return df["charttime_offset"].astype("float64")
    if "offset_min" in df: return df["offset_min"].astype("float64")
    raise KeyError(df.columns.tolist())

def load(db, var, stays=None, extra=()):
    f = f"{P}/{db}_{var}.parquet"
    cols = [c for c in pq.ParquetFile(f).schema_arrow.names if c in ("stay_id", "charttime", "charttime_offset", "offset_min", var, *extra)]
    filt = ds.field("stay_id").isin(list(stays)) if stays is not None else None
    d = ds.dataset(f).to_table(columns=cols, filter=filt).to_pandas()
    d["t"] = tmin(d, db); d = d.drop(columns=[c for c in ("charttime", "charttime_offset", "offset_min") if c in d])
    if var != "ventilator_mode_raw": d[var] = pd.to_numeric(d[var], errors="coerce").astype("float32")
    d["stay_id"] = d["stay_id"].astype("int64")
    return d.dropna(subset=["t"]).sort_values("t")

def asof(anchor, d, var, name=None):
    name = name or var
    d = d[["stay_id", "t", var]].dropna().rename(columns={var: name, "t": f"t_{name}"})
    m = pd.merge_asof(anchor.sort_values("t"), d.sort_values(f"t_{name}"), left_on="t", right_on=f"t_{name}",
                      by="stay_id", direction="backward")
    m[f"lag_{name}"] = (m["t"] - m[f"t_{name}"]).astype("float32")
    return m.drop(columns=[f"t_{name}"])

def cohort(db):
    c = pd.read_parquet(f"{P}/cohort_{db}.parquet")
    age = pd.to_numeric(c["anchor_age"].astype(str).str.replace(">", "").str.strip(), errors="coerce")
    sexcol = "gender" if db == "mimic" else ("sex" if "sex" in c and c["sex"].notna().any() else "gender")
    sex = c[sexcol].astype(str).str.upper().str[0].where(lambda s: s.isin(["M", "F"]))
    c = c.assign(age=age, sex=sex)
    n0 = len(c)
    c = c[c.age >= 18]
    n_adult = len(c)
    # first ICU stay per patient
    if db == "eicu":
        c = c[c.unitvisitnumber == 1].sort_values("stay_id").drop_duplicates("subject_id")
    else:
        c = c.sort_values("icu_intime").drop_duplicates("subject_id")
    # height -> PBW (first recorded; HiRID imputed heights excluded)
    if db == "mimic":
        h = pd.read_parquet(f"{P}/mimic_height_cm.parquet", columns=["stay_id", "charttime", "height_cm"]).sort_values("charttime")
        h = h.drop_duplicates("stay_id")[["stay_id", "height_cm"]]
    else:
        cols = ["stay_id", "admission_height_cm"] + (["is_imputed"] if db == "hirid" else [])
        h = pd.read_parquet(f"{P}/{db}_admission_height_cm.parquet", columns=cols)
        if db == "hirid": h = h[~h.is_imputed.astype(bool)]
        h = h.rename(columns={"admission_height_cm": "height_cm"}).drop_duplicates("stay_id")[["stay_id", "height_cm"]]
    h["height_cm"] = pd.to_numeric(h.height_cm, errors="coerce")
    c = c.merge(h, on="stay_id", how="left")
    c.loc[~c.height_cm.between(120, 220), "height_cm"] = np.nan
    c["pbw"] = np.where(c.sex == "M", 50, 45.5) + 0.91 * (c.height_cm - 152.4)
    c.loc[c.sex.isna() | c.height_cm.isna(), "pbw"] = np.nan
    c["stay_id"] = c.stay_id.astype("int64")
    return c[["stay_id", "subject_id", "age", "sex", "height_cm", "pbw"]], dict(n_cohort_file=n0, n_adult=n_adult, n_first_stay=len(c))

def build(db):
    c, att = cohort(db)
    stays = set(c.stay_id)
    pl = load(db, "plateau_pressure", stays)
    att["plateau_rows"] = len(pl); att["plateau_stays"] = pl.stay_id.nunique()
    ev = pl[["stay_id", "t", "plateau_pressure"]]
    for var, nm in [("peep_total", "peep"), ("tidal_volume_observed", "vt_obs"), ("tidal_volume_set", "vt_set"),
                    ("respiratory_rate_set_vent", "rr_set")]:
        ev = asof(ev, load(db, var, stays), var, nm)
    if db == "hirid":
        # SUPER-ICU's HiRID mode labels do not match HiRID ordinal_vars_ref.csv (itemid 3845); map raw codes officially.
        mo = load(db, "ventilator_mode_raw", stays, extra=("mode_code",))
        hm = pd.read_csv(os.environ.get("HIRID_MODE_MAP", "hirid_mode_official_map.csv"))
        mo["mode_class"] = mo.mode_code.map(dict(zip(hm.mode_code, hm.mode_class))).fillna("unknown")
    else:
        mo = load(db, "ventilator_mode_raw", stays)
        mo["mode_class"] = mo.ventilator_mode_raw.astype(str).map(dict(zip(MODE_MAP.mode_raw.astype(str), MODE_MAP.mode_class))).fillna("unknown")
    ev = asof(ev, mo, "mode_class")
    rr = load(db, "respiratory_rate", set(ev.stay_id))
    ev = asof(ev, rr, "respiratory_rate", "rr_tot")
    ev = ev.merge(c[["stay_id", "pbw"]], on="stay_id", how="left")
    ev.to_parquet(f"{OUT}/events_{db}.parquet", index=False)
    c.to_parquet(f"{OUT}/cohort_{db}_pbw.parquet", index=False)
    att.update(db=db, event_rows=len(ev), event_stays=ev.stay_id.nunique())
    return att

if __name__ == "__main__":
    rows = []
    for db in DBS:
        rows.append(build(db)); print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/attrition_events.csv", index=False)
