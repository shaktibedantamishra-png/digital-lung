# Post hoc analysis addendum v1 (NOT pre-specified; written 5 Oct 2026 before any result below was computed)

Trigger: external presubmission audit (findings 1, 2, 3, 8 and the hysteresis item). Locked spec analysis_spec_v1.md is unchanged; every analysis here is additional and is labelled "not pre-specified" in the manuscript. Base data: primary filtered events (tolerance 60 min, exhaled Vt, permissive passivity), phenotype definitions exactly as S6 unless stated. Known before writing (descriptive, not outcome): PEEP and Vt share the plateau timestamp in 98.3% (MIMIC-IV), 99.5% (eICU-CRD) and 20.8% (HiRID; 90th percentile lag 1.9 min) of primary events. Seed 20261005.

## A1 Join-lag artefact (finding 1)
- A1a Lag strata: per event, max(lag_peep, lag_vt) in {0, >0–5, >5–15, >15–60 min}. Phenotype recomputed using only events in each stratum (stratum must independently meet S6 eligibility).
- A1b Same-time sensitivity: MIMIC-IV and eICU-CRD: lag_peep = 0 and lag_vt = 0; HiRID: both ≤5 min.
- A1c Ambiguity exclusion: drop events where the next charted PEEP after the plateau time (within 60 min) differs from the as-of PEEP by ≥1 cmH2O.

## A2 Measurement-noise null calibration (finding 1)
Per phenotype-eligible stay and qualifying Vt bin, keep observed event times, PEEP values and Vt values. Null plateau: Pplat* = PEEP_true + Vt / C_stay + e, where C_stay = stay median Crs over primary events, e ~ N(0, s_e^2), s_e = pooled within-stay SD of Pplat among events sharing stay, Vt bin and rounded PEEP level (replicate SD; pooled per database). PEEP noise: PEEP_true = PEEP_obs − u, u ~ N(0, s_u^2); s_u per database = pooled within-stay SD of (PEEP total − PEEP set), as-of PEEP set within 60 min (MIMIC-IV, HiRID); eICU-CRD uses the larger of the two. Null slope is 1 in expectation except for observed Vt covariation. 200 replicates. Outputs: null share <0.9 and >1.1 (mean, 2.5–97.5 percentile); excess = observed share − null mean. Variants: s_u = 0 (noise in Pplat only) and s_u doubled.
- A2b Set-PEEP predictor (MIMIC-IV, HiRID): phenotype with as-of PEEP set (≤60 min) replacing PEEP total; set PEEP has no measurement error by definition.

## A3 Mode (finding 2)
Stay-level mode stratum = modal raw mode across the stay's primary events (mode within 720 min, D3).
- MIMIC-IV: volume control constant flow = CMV/ASSIST, CMV, VOL/AC, (S) CMV; volume-targeted pressure-regulated = CMV/ASSIST/AutoFlow, CMV/AutoFlow, PRVC/AC, APV (cmv); pressure control = PCV+Assist, PCV+, PRES/AC, P-CMV; other/unknown = remainder.
- HiRID (official codes): CMV = volume control; DGV = pressure control; Druckreg. Vol.Kontr. = volume-targeted pressure-regulated; other/unknown = remainder. eICU-CRD: not applicable (D4).
- A3b Covariate-adjusted slope: per stay, OLS Pplat ~ PEEP + Vt/PBW over all primary events with PBW (no binning); eligibility ≥2 rounded PEEP levels, span ≥3 cmH2O, ≥5 events, Vt not collinear (|r(PEEP, Vt)| < 0.95); slope = PEEP coefficient, classified with S6 thresholds.
- A3c Hysteresis restricted to steps whose both events have mode class 'controlled' (MIMIC-IV, HiRID).

## A4 Population subgroups (finding 3)
- P/F at the stay's first qualifying event (t_first): PaO2 nearest within ±4 h; FiO2 as-of ≤4 h before the PaO2 (fraction; values >1 divided by 100; 0.21–1.0 kept); PEEP at t_first ≥5. Strata: P/F ≤300 (hypoxaemic) and ≤150 (moderate–severe); >300; missing.
- Care unit: MIMIC-IV surgical = CVICU, SICU, TSICU, Neuro SICU; medical = MICU, CCU; mixed = MICU/SICU; other = remainder. eICU-CRD surgical = SICU, CSICU, CTICU; medical = MICU, Cardiac ICU; mixed = Med-Surg ICU, CCU-CTICU; other = Neuro ICU. HiRID: diagnosis not available in the harmonised data; not analysed.
- Time window stratum already reported (Table S4).

## A5 eICU-CRD between-hospital heterogeneity (finding 8)
Hospitals with ≥20 phenotype stays. Logit proportion recruitment-consistent; DerSimonian–Laird random effects; pooled proportion with 95% CI, tau^2, I^2, 95% prediction interval (t, k−2 df). Also the overdistension-consistent share.

## Reporting rule
All results reported whatever their direction. If the excess over the null (A2) is small, the phenotype claim in the manuscript is weakened accordingly.
