# Analysis specification v1 — locked 2026-09-25 before any outcome computation
Databases: MIMIC-IV, eICU-CRD, HiRID (SUPER-ICU harmonised parquets). AmsterdamUMCdb excluded (plateau variable maps to a non-measured channel; S1).

## S1 Cohort
Adults (age >= 18 y), first ICU stay per patient where identifiable, with >= 1 plateau pressure record. Source: cohort_<db>.parquet.

## S2 Event construction (plateau-anchored)
Anchor = each plateau pressure record. Backward as-of join by stay (most recent value at or before the plateau time) of: PEEP total, Vt observed (exhaled; primary) and Vt set (sensitivity), ventilator mode, RR total, RR set.
Tolerance: 60 min primary; 120 and 240 min sensitivity. Values older than the tolerance -> missing.
Plausibility: Pplat 5-50 cmH2O; PEEP 0-30 cmH2O; Vt 150-1200 mL; Pplat > PEEP.

## S3 Passivity
Primary (permissive): exclude events whose most recent mode is classified spontaneous (CPAP/PSV/SPONT etc.; mapping table mode_mapping.csv built from raw strings before outcomes). Unknown mode retained.
Sensitivity (strict): controlled modes only. Sensitivity (RR rule): RR total <= RR set + 2 breaths/min where both available.
eICU: if mode coverage is < 10% of events, only the permissive definition is applicable; stated as a limitation.

## S4 PBW
Devine: male 50 + 0.91(height_cm - 152.4); female 45.5 + 0.91(height_cm - 152.4). Height 120-220 cm; first recorded height. Stays without sex or height excluded from PBW-dependent analyses only.

## S5 Static compliance
Crs = Vt / (Pplat - PEEP) [mL/cmH2O] per event (Vt observed, primary). Per-stay median. Report per database n, median, IQR.

## S6 PEEP-response phenotype
Within stay, events grouped into Vt bins: floor(Vt / PBW) in 1 mL/kg steps (primary); 50 mL absolute bins (sensitivity).
A bin qualifies if it contains >= 2 distinct PEEP levels spanning >= 3 cmH2O. Bin slope = OLS slope of Pplat on PEEP. Stay slope = median of qualifying-bin slopes.
Classes: recruitment-consistent < 0.9; neutral 0.9-1.1; overdistension-consistent > 1.1.
Time window (descriptive): hours from first plateau record to the first event in the qualifying bin: <2, 2-6, 6-24, >24 h.

## S7 Hysteresis (direction asymmetry)
Chronologically consecutive events within stay, same Vt bin, |dPEEP| >= 2 cmH2O, gap <= 4 h. Step slope = dPplat/dPEEP.
Per stay: median slope for increases and for decreases. Paired stays = both directions present. Effect = median within-stay difference (decrease - increase); Wilcoxon signed-rank p; 2000-replicate percentile bootstrap CI (resampling stays).

## S8 Statistics
Proportions: Wilson 95% CI. Medians: 2000-replicate percentile bootstrap. Between-database phenotype distribution: chi-square, Cramer's V with bootstrap CI. Seed 20260925. No equivalence claims.

## S9 Reporting rule
Results are reported as computed. Definitions may be changed only for data-structure reasons found in the audit (units, field availability), documented in the deviations log below, never in response to outcome values.

## Deviations log (all data-structure driven; logged before outcome computation)
D1 Cohort files are SUPER-ICU pre-filtered: all rows adult; eICU and HiRID already one stay per patient. MIMIC: first ICU stay per subject by icu_intime.
D2 HiRID ventilator mode: SUPER-ICU labels for itemid 3845 do not match HiRID ordinal_vars_ref.csv (e.g. code 3 = CMV labelled 'simv'; code 8 = DU/ST labelled 'volume_control'; code 1 = Stand By labelled 'spontaneous'). Raw codes re-mapped from the official reference (hirid_mode_official_map.csv). Under the official map 99.2% of HiRID plateau rows with a mode are CMV.
D3 Mode as-of tolerance 720 min (not 60): mode is a setting charted intermittently (MIMIC only 52% of plateau rows have a mode within 12 h). PEEP, Vt, RR keep the 60 min primary tolerance.
D4 eICU: ventilator-mode table covers <1% of plateau rows (0 within 12 h). Permissive passivity is therefore no filter in eICU; strict policy not applicable; RR rule is the only passivity sensitivity for eICU.
D5 HiRID imputed heights (is_imputed) excluded from PBW.
D6 PEEP levels for the >=2 levels / >=3 cmH2O span rule are counted on PEEP rounded to 1 cmH2O (HiRID PEEP total is a measured continuous value; 85% non-integer). OLS slope uses unrounded values.
