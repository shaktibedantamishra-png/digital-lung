# Digital-lung trial engine: specification v1 (DRAFT, for sign-off before any run)

Status: nothing has been coded or run under this spec. It replaces the original engine, whose source code could not be located. It was written **with knowledge of the original engine's outputs**, and the manuscript must disclose this. No parameter may be changed after the first run except for a coding error, which is logged with its reason.

Every parameter is labelled as one of:
- **[data]**: calibrated from the locked pipeline;
- **[lit-abstract]**: value taken from the cited abstract, verified this session;
- **[assumption]**: no data constraint; it gets a sensitivity range.

## E1. Virtual population
- Source: per-stay C_rs (compliance_perstay.parquet) joined to PBW (pipe_out/cohort_<db>_pbw.parquet), with complete cases only. **[data]**
- Primary: 2000 (C_rs, PBW) pairs per database are resampled with replacement from the observed pairs. This is empirical resampling, so there is no copula and no tail attenuation. Seed 20260929.
- Sensitivity: a Gaussian copula on ranks with empirical marginals.
- Uncertainty: 200 bootstrap re-draws of the population. Report how often each minimum-cost setting is selected.

## E2. Lung model (one non-recruitable and one recruitable compartment)
- Aerated baby-lung volume at ZEEP: FRC = EL_spec × C_L, where C_L = C_rs / κ.
  - EL_spec = 13.5 cmH2O **[lit-abstract: Chiumello 2008, ARDS group 13.5 ± 4.1]**.
  - κ = E_L/E_rs = 0.7 **[assumption; abstract range 0.33–0.95]**.
- Recruitable gas volume: V_R = r × FRC, with r = 0.3 **[assumption; sensitivity 0.1, 0.5]**.
- Each recruitable unit has a closing pressure TCP ~ N(μ_c, σ²) and an opening pressure TOP = TCP + g.
  - μ_c = 5 **[lit-abstract: Crotti 2001, closing-pressure mode 5]**.
  - g = 15, so the TOP mode is 20 **[lit-abstract: Crotti 2001, opening-pressure mode 20]**.
  - σ = 5 **[assumption; sensitivity 3, 8]**. The abstract states that the distributions are normal but gives no SD.
- In steady state at a given PEEP and P_plat (= PEEP + V_T/C_rs), each unit is in one of three states:
  - open throughout if TCP < PEEP and TOP ≤ P_plat;
  - cyclic (collapses at end-expiration, reopens each breath) if TCP ≥ PEEP and TOP ≤ P_plat;
  - closed throughout if TOP > P_plat.
- f_open and f_cyc are the corresponding normal-probability fractions, computed in closed form.

**Decision needed (D-A): does compliance depend on PEEP?**
- **A (recommended for this paper):** MP uses the measured C_rs, fixed. Recruitment state enters only the penalties. MP therefore stays a purely data-calibrated quantity. The engine then predicts ΔP_plat/ΔPEEP = 1 and cannot reproduce the observed phenotype; we state this.
- **B:** C_rs(PEEP) scales with the aerated volume, (FRC + V_R·f_open), calibrated at each stay's median clinical PEEP. The engine can then predict the phenotype slope distribution, giving a real output-validation check. The cost is extra assumptions, and the phenotype would partly calibrate r and σ, which risks overclaiming recruitability. Better suited to the validation project.

## E3. Mechanical power **[identical to the planner module and to the stored outputs]**
MP = 0.098 · RR · [0.5 · E_rs · V_T² + RR · K · R_aw · V_T² + PEEP · V_T]
- E_rs = 1/C_rs (cmH2O/L); K = (1 + I:E)/(60 · I:E) with I:E = 0.5.
- R_aw = 10 cmH2O·s/L **[assumption; sensitivity 5, 15]**.
- V_T in L. MP in J/min.

## E4. Penalties (J/min; energy per breath × RR, so the scaling with rate is physical and explicit)
- **Cyclic recruitment/derecruitment:** P_dr = 0.098 · RR · V_R · f_cyc · (μ_c + g). This is the energy to reopen cyclically collapsing units at their modal opening pressure.
- **Overdistension:**
  - Total strain ε = (C_rs·PEEP + V_T) / (FRC + V_R·f_open). V* = ε*·(FRC + V_R·f_open).
  - P_od = 0.098 · RR · 0.5 · E_L · [V_ei² − max(V_PEEP, V*)²], counted only while V_ei > V*. Here V_PEEP = C_rs·PEEP, V_ei = V_PEEP + V_T and E_L = κ·E_rs. This is the tidal elastic energy delivered above the strain threshold.
  - ε* = 1.5 **[lit-abstract: Protti 2011, pigs, damage above 1.5–2; sensitivity 2.0]**.
- Composite cost: J = MP + w_dr·P_dr + w_od·P_od, with w ∈ {0, 0.5, …, 4}² (81 combinations). J is minimised over the population mean.
- Required properties, tested in E7: P_od is non-decreasing in PEEP and V_T; P_dr is ≥ 0 and exactly proportional to RR.

## E5. Strategy grid and analyses
1. **Main grid:** V_T 4, 6, 8, 10 mL/kg PBW × PEEP 5, 8, 10, 12, 15, 18, 20 cmH2O × RR 15, 20, 25/min (84 strategies). The ranked-MP panel uses all 84 strategies.
2. **Matched minute ventilation:** V_E = 100 mL/kg PBW/min, RR = V_E/V_T (limited to 10–45), PEEP 10.
3. **Matched alveolar ventilation (new; this is the metabolic constraint):** V_A = RR·(V_T − V_D) held at the value for 6 mL/kg at 20/min, with V_D = 2.2 mL/kg PBW **[assumption; sensitivity 1.5, 3.0]**.
4. **Protective-limit variant:** exclude strategies with P_plat > 30 or ΔP > 15 cmH2O.
5. **One-at-a-time sensitivity:** r, σ, κ, ε*, R_aw and V_D. Report the minimum-cost V_T, PEEP and RR in each run.

## E6. Reporting rules
- Minimum-cost settings are minima of a surrogate within a closed grid, never recommendations.
- Report every weight combination and every sensitivity run, not a selected subset.
- If a result contradicts the v12 text, the text changes. The spec does not.

## E7. Verification (unit tests; all must pass before any trial run)
1. The MP function matches a hand calculation for 3 cases to 1e-9.
2. With w = 0, the matched-MV MP reproduces the closed-form expression.
3. f_open + f_cyc + f_closed = 1 on a grid.
4. P_od is monotonic non-decreasing in PEEP and V_T on a 0.5-step fine grid.
5. P_dr(RR) / RR is constant.
6. Seeds are reproducible: identical outputs on rerun.

## E8. Deliverables
- engine/digital_lung.py and tests/test_engine.py, with SHA-256 of the spec and code logged at the first run.
- trial_main.csv, trial_matched_mv.csv, trial_matched_va.csv, weight_sweep.csv and sensitivity_oat.csv.
- New Figure 1 panels B–D and a rewritten S1, with every equation copied from code.
- Manuscript changes: in-silico Results, Discussion correction paragraph, Table 1 (code verification = Done) and a methods disclosure.

## References for parameters (DOIs checked on CrossRef this session)
- Chiumello D et al. Am J Respir Crit Care Med 2008;178:346–55. doi:10.1164/rccm.200710-1589OC
- Protti A et al. Am J Respir Crit Care Med 2011;183:1354–62. doi:10.1164/rccm.201010-1757OC
- Crotti S et al. Am J Respir Crit Care Med 2001;164:131–40. doi:10.1164/ajrccm.164.1.2007011
- Hickling KG. Am J Respir Crit Care Med 1998;158:194–202. doi:10.1164/ajrccm.158.1.9708049 (compartment-threshold modelling precedent)
