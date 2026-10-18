# digital-lung-icmx

Code, specifications and aggregate results for "What routine ICU mechanics data can and cannot constrain in a human-data-bounded digital lung" (Intensive Care Medicine Experimental, submission candidate).

No patient-level data are included. The source databases (MIMIC-IV, eICU-CRD, HiRID) are available to credentialed users via PhysioNet. Pipeline scripts read harmonised SUPER-ICU parquet files that are not redistributed; the path is set through the environment variable SUPERICU_DATA.

## Versions and dates

| Document | Version | Date | Status |
|---|---|---|---|
| spec/analysis_spec_v1.md | v1 | 25 Sep 2026 (file date) | Locked analysis plan for the primary analyses |
| spec/posthoc_addendum_v1.md | v1 | 5 Oct 2026 | Written before the additional analyses were computed; SHA-256 in spec/posthoc_addendum_v1.sha256 |
| spec/posthoc_deviations.md | log | 7–18 Oct 2026 | Deviations PD1–PD8, logged as they occurred |
| spec/engine_spec_v1.md | v1 | not dated in file | Header reads "DRAFT, for sign-off before any run"; see note below |

Note on the engine specification: the manuscript describes the engine specification as locked before the engine was run. The header of spec/engine_spec_v1.md still says DRAFT, and no lock hash file is included. Either add the lock record (hash and date) or change the manuscript wording before release.

## Contents

- spec/ — analysis specification, engine specification, post hoc addendum and deviation log
- pipeline/ — event construction (build_events.py), filters, phenotype classification, hysteresis steps
- engine/ — virtual-population trial engine (digital_lung.py), trial runner (run_trials.py), joint Latin-hypercube sensitivity (lhs_sweep.py)
- tests/ — unit tests for the engine (six pre-specified tests)
- analysis/ — null simulation, excess bootstrap CIs, PD6/PD8 hysteresis (null_bootstrap.py); compliance with set PEEP (compliance_setpeep.py); shared setup (common.py)
- results/ — aggregate tables only (no stay-level rows); hysteresis_setpeep_pd6_pd8.csv holds PD6 and PD8 hysteresis
- requirements.txt — package versions used

## Cohort extraction (external)

The SUPER-ICU harmonisation step that produces the per-database parquet files is NOT included in this repository. It is an external pipeline; the repository starts from its output. Pipeline inputs: harmonised SUPER-ICU parquet files for MIMIC-IV (v3.1), eICU-CRD (v2.0) and HiRID (v1.1.1), accessed via PhysioNet. Consequently the repository cannot regenerate the parquet inputs from raw PhysioNet downloads.

## Reproduction

1. Set SUPERICU_DATA to the harmonised parquet directory (note: pipeline/build_events.py reads this variable; the path in the original build was redacted for release).
2. python pipeline/build_events.py
3. Run the phenotype and hysteresis functions from pipeline/ on the events.
4. python tests/test_engine.py
5. python engine/run_trials.py ; python engine/lhs_sweep.py
6. python analysis/null_bootstrap.py (noise-only null, PEEP-error variant, excess CIs, PD6/PD8 hysteresis)
7. python analysis/compliance_setpeep.py

Environment: see requirements.txt. Scripts 6 and 7 write to results_reproduced/ when run from the repository root; the copies in results/ were produced by this procedure and match it exactly.

Analysis parameters: seed 20260929 (engine), 20261005 (post hoc), 20261006 and 20261007 (null simulations and bootstraps).

## Known limitations of this release

- The engine is an exploratory surrogate. Its minimum-cost settings are not clinical recommendations.
- eICU-CRD records set PEEP only; intrinsic PEEP cannot be assessed there.
- The no-intrinsic-PEEP subsets (MIMIC-IV n=593, HiRID n=286) are small and selected; hysteresis in these subsets is exploratory.
- Author and ethics placeholders in the manuscript are not part of this repository.

## Licence and citation

Code: to be chosen by the authors before release (not yet applied). Cite the Zenodo DOI for the release once minted.
