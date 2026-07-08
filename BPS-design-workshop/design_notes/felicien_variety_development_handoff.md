# Felicien Variety Development Scheme Handoff

Design-only handoff for later BPS 0.2.0 implementation. Updated against `Felicien_Variety_development_v2.jpeg`.

No R code was written, no R commands were run, and no validation or smoke tests were performed during this session. This document records the approved biological interpretation, event/scheduler contract, and implementation plan for a future coding session.

## Status

- Mode: design-only
- Skill used: breeding-scheme-drafter
- Gate 1 interpretation: approved
- Gate 2 event verbs and scheduler contract: approved
- Gate 3 implementation plan: drafted for handoff
- v2 image review: incorporated before implementation
- Implementation: deferred
- Validation: deferred

## Biological Objective

This is a fast-exploitation hybrid rice variety-development scheme focused on one elite R-line family. The recurrent pipeline improves one R-line pool, while fixed female testers represent the opposite heterotic group. The released product is hybrids only.

The scheme intentionally exploits high within-family relatedness among elite R lines selected from one outstanding biparental family. The goal is faster genomic selection and hybrid release, not long-term genetic diversity maintenance.

## Approved Gate 1 Interpretation

### Starting Material

- Start from 5 fixed elite R-line founders.
- All 5 founders were selected from the same previous biparental family.
- The historical biparental selection does not need to be simulated in the first implementation; record it as provenance unless a later experiment studies that history.

### Training Phase (`Cycle 1: Training` in v2 image)

1. Use the 5 fixed elite R-line founders.
2. Make 10 elite x elite crosses among them by default. Keep the implementation configurable if later experiments need 5-10 crosses, but the v2 diagram labels this count as 10.
3. Generate 800 total R-line descendants across crosses by default. Keep the implementation configurable if later experiments need 800-1000, but the v2 diagram labels F2, F3, F4, F5, F6, and MET testcross scale as 800.
4. Advance the same cohort by SSD/RGA with no selection from F1 through F5. F2-F4 are temporary RGA intermediates.
5. At F5, continue RGA/seed increase and genotype the training R lines; the v2 image marks F5 as the point for `G` and `GS`.
6. Advance the training R-line cohort to F6 for testcross production.
7. Use fixed female testers. The earlier interpretation used 2 testers; v2 still shows a tester block but does not specify tester count in training, so keep tester count configurable and default to 2 for consistency with final-cycle arithmetic.
8. Produce training testcross hybrids for MET. The v2 image labels the MET testcross stage as 800, so the default simulation target is 800 connected testcross entries/plots instead of the earlier capped sparse-MET assumption.
9. Run partner MET on testcrosses for R lines and testers.
10. Estimate R-line GCA values across traits, with SCA retained as reporting/diagnostic information if the implemented design supports it.
11. Construct a desired-gain selection index from GCA values.
12. Select the best 5 R lines to start recurrent GS.
13. Train one fixed GP model from the training-cycle genotype and GCA/index information.

### Recurrent GS Phase

There are five recurrent GS cycles.

Cycles 1-4:

1. Use the current selected 5 R lines directly as parents.
2. Make 10 elite x elite crosses by default, matching the v2 label. Keep this configurable.
3. Generate 50 total candidates across all crosses.
4. Advance by SSD/RGA to F5 with no field phenotyping; F2-F4 are temporary RGA intermediates.
5. Genotype all 50 F5 candidates.
6. Predict index GEBV using the unchanged GP model.
7. Select the best 5 F5 candidates.
8. Use those 5 immediately as parents for the next cycle; the v2 image explicitly marks this as a 4-cycle recurrent loop.

Cycle 5:

1. Make recurrent crosses from the cycle 4 parents, with 10 elite x elite crosses by default.
2. Generate and advance 50 total candidates to F5.
3. Genotype all 50 F5 candidates.
4. Predict index GEBV using the fixed model.
5. Select 30 immediately at F5.
6. Advance only those 30 lines to F6.
7. Cross the 30 selected R lines to fixed female testers. With 2 testers, this produces 60 AYT hybrid combinations.
8. Send the 60 hybrids to partner AYT.
9. Release all hybrids whose observed AYT yield strictly outperforms the local commercial check.

### Trait Map

Training GCA/index traits:

- yield
- blast
- hoja blanca
- Bulkhoderia
- milling yield
- milling quality
- white center

Recurrent GS prediction target:

- desired-gain/index GEBV trained from training-cycle GCA/index values

Final AYT release criterion:

- observed hybrid yield strictly outperforming local commercial check

Other AYT traits:

- measured and reported
- not release-gating in the first design

### Agreed Simplifications

- Drop standalone F6 disease screening.
- Do not model F6/testcross/AYT in recurrent cycles 1-4.
- Do not retrain the GP model during the five recurrent cycles.
- Keep F2-F4 as temporary SSD/RGA intermediates.
- Focus recurrent parent selection on GCA/index, not SCA.
- Store AYT data for possible future model updating in a later round, but do not use it for retraining in this five-cycle run.

## Persistent And Temporary Populations

### Store In BPS State

- `R_Founders`: 5 fixed elite R-line founders.
- `Training_RLines_F5`: 800 genotyped training R lines at F5 by default.
- `Training_Testcross_MET`: connected R-line x tester hybrids evaluated in partner MET.
- `Training_RLines_GCA`: training R lines carrying GCA, SCA, and desired-gain index summaries.
- `GS_Model`: fixed GP model trained once from training GCA/index.
- `R_Parents_Cycle_0`: best 5 R lines from training cycle.
- `R_Candidates_F5_Cycle_i`: 50 genotyped F5 candidates for recurrent cycle `i`.
- `R_Parents_Cycle_i`: selected 5 parents for recurrent cycles 1-4.
- `R_Final30_F5`: 30 selected F5 lines from cycle 5.
- `R_Final30_F6`: final 30 R lines after continued SSD/RGA/seed increase to F6.
- `AYT_Hybrids`: 60 hybrids from 30 R lines x 2 testers.
- `Released_Hybrids`: hybrids outperforming the local commercial check for yield.

### Keep Temporary Unless Later Needed

- F1/F2/F3/F4 intermediate generations.
- Unselected recurrent F5 candidates after reporting.
- Testcross production intermediates before MET/AYT entries.
- Plot-level allocation objects, unless implementation needs explicit plot records.

## Gate 2 Event Verbs

### `initialize_R_founders(state, cfg, year)`

Seeds the 5 fixed elite R-line founders. No selection occurs inside the scheme unless a later experiment simulates the historical biparental family.

### `make_training_R_crosses(state, cfg, year)`

Uses `R_Founders` to make configurable 5-10 elite x elite crosses. Produces the training family base.

### `advance_training_RLines_by_SSD(state, cfg, year)`

Advances training progeny by SSD/RGA with no selection. Stores the F5 training cohort as `Training_RLines_F5`; the v2 image then carries this cohort forward for F6 testcross production.

### `genotype_training_RLines(state, cfg, year)`

Genotypes all training R lines.

### `make_connected_training_testcrosses(state, cfg, year)`

Creates connected R-line x tester hybrid entries using fixed testers. The v2 image labels the training MET testcross scale as 800; tester count remains configurable, with 2 as the current default pending confirmation.

### `run_partner_MET_and_estimate_GCA(state, cfg, year)`

Runs partner MET on the connected training testcross set. The v2 image no longer shows the earlier 200-plots-per-location cap, so default implementation should target the diagram scale of 800 MET testcross observations/entries unless a later cfg explicitly requests per-location caps. Stores hybrid MET entries and attaches GCA, SCA, and index summaries back to R lines.

### `select_training_R_parents_by_GCA_index(state, cfg, year)`

Selects best 5 R lines using desired-gain index on GCA across all traits. SCA remains reportable/diagnostic.

### `train_fixed_GS_model_from_GCA_index(state, cfg, year)`

Trains one GP model from genotypes plus GCA/index response. This model stays unchanged during the five recurrent cycles.

### `make_recurrent_R_crosses(state, cfg, year, cycle)`

Uses current 5 selected R parents to make 10 crosses by default for the requested cycle; keep the number configurable.

### `advance_recurrent_to_F5_and_genotype(state, cfg, year, cycle)`

Generates 50 total candidates across all crosses, advances by SSD/RGA to F5, and genotypes all 50.

### `predict_and_select_recurrent_R_parents(state, cfg, year, cycle)`

For cycles 1-4, predicts index GEBV using the fixed GP model and selects the best 5 F5 candidates as next-cycle parents.

### `predict_and_select_final30_RLines(state, cfg, year, cycle = 5)`

For cycle 5, predicts index GEBV at F5 and selects 30 immediately.

### `advance_final30_F5_to_F6(state, cfg, year)`

Continues only the 30 selected lines from F5 to F6.

### `make_final_testcross_hybrids(state, cfg, year)`

Crosses all 30 selected R lines to the same 2 fixed testers, producing 60 candidate hybrid combinations.

### `run_partner_AYT(state, cfg, year)`

Runs partner AYT on the 60 hybrids for the same trait set. AYT data are stored for future model updating, but do not update the model in this five-cycle run.

### `release_hybrids_from_AYT(state, cfg, year)`

Releases every hybrid whose observed AYT yield strictly outperforms the local commercial check. Other traits are reported but do not gate release in this version.

### `record_yearly_outputs(results, results_base, state, year, cfg)`

Reports stage sizes, mean genetic values, index values, prediction accuracy where available, GCA/SCA summaries, MET/AYT trait summaries, released-hybrid count, and costs.

## Approved Scheduler Contract

1. Initialize 5 R founders.
2. Build training crosses.
3. Advance 800 training R lines by SSD/RGA to F5 by default.
4. Genotype training R lines at F5 and continue them to F6 for testcross production.
5. Produce connected training testcrosses.
6. Run partner MET at the v2 default testcross scale of 800 and estimate GCA, SCA, and GxE where supported.
7. Select 5 R parents by GCA desired-gain index.
8. Train one fixed GS model.
9. Run recurrent cycles 1-4:
   - make recurrent crosses from current 5 parents
   - advance 50 candidates to F5
   - genotype all 50
   - predict index GEBV
   - select 5 as next parents
   - record cycle results
10. Run recurrent cycle 5:
    - make recurrent crosses from cycle 4 parents
    - advance 50 candidates to F5
    - genotype all 50
    - predict index GEBV
    - select 30
    - advance selected 30 to F6
    - cross 30 R lines to 2 testers
    - run partner AYT on 60 hybrids
    - release hybrids outperforming local yield check
11. Store AYT data for possible future model update, but do not retrain within this run.

## Gate 3 Implementation Plan

### Verified BPS/AlphaSimR Patterns

Use these with high confidence later:

- Store cohorts with `put_stage_pop()`.
- Update stored cohorts with `bp_update_stage_pop()`.
- Select available cohorts with `select_latest_available()`.
- Skip expected warm-up gaps with `bp_skip_if_no_input()`.
- Make R-line crosses with AlphaSimR `randCross()` or `makeCross()`.
- Advance SSD/RGA with repeated AlphaSimR `self()`, keeping F2-F4 temporary.
- Genotype stored cohorts with `run_genotyping()`.
- Run multi-location trials with `run_phenotype_trial()`, including `use_env_control`, `env_means`, `env_mean_sd`, `env_year_sd`, and per-environment logs.
- Use BPS synthetic trait/index machinery for desired-gain index selection where it fits.
- Train and predict GP using BPS GP patterns around `run_train_gp_model()`, `predict_ebv_pop()`, or explicit `AlphaSimR::RRBLUP()` model storage.

### Custom Or To-Verify Later

These are biologically central but need careful implementation checking:

- Connected training testcross design across fixed testers, defaulting to the v2 scale of 800.
- MET allocation for the v2 diagram scale of 800; only apply per-location caps if later requested by cfg.
- Mixed-model summaries for GCA, SCA, GxE, and environment correlations.
- Attaching GCA/SCA/index summaries back onto R-line `pop@misc`.
- Training a multi-trait or index-response GP model from GCA-derived values.
- Final release rule comparing each AYT hybrid to local commercial checks.
- Confirm whether the training tester count should remain 2; v2 only labels the final AYT arithmetic (`30 x 2 = 60`) explicitly.

### Event Implementation Notes

`initialize_R_founders`

- Seed or select 5 elite R-line founders from the initialized simulation state.
- Store as `R_Founders`.
- Record same-family origin as provenance.

`make_training_R_crosses`

- Use `R_Founders` to make `cfg$n_training_crosses` among the 5 parents.
- Cross plan should be explicit and configurable: half-diallel-like, random among founders, or user-specified.
- Store only if needed; otherwise feed directly into SSD.

`advance_training_RLines_by_SSD`

- Generate `cfg$n_training_rlines` total across crosses.
- Default v2 production scale: 800. Earlier 800-1000 range can remain as an experiment parameter.
- Repeated selfing advances to the inbred stage.
- No selection.
- Store `Training_RLines_F5`.
- Log line-development and seed-increase costs.

`genotype_training_RLines`

- Run genotyping on all `Training_RLines_F5`.
- Required cfg fields: chip, cost, duration.
- Cohort persists because GP training and GCA summaries depend on it.

`make_connected_training_testcrosses`

- Create hybrid entries from training R lines x fixed testers.
- Default v2 MET testcross scale is 800.
- Prefer a connected design; later code may physically create only assigned hybrids or create full combinations and subset to the design.
- Store `Training_Testcross_MET`.

`run_partner_MET_and_estimate_GCA`

- Run partner MET with the approved trait list.
- Implement MET allocation separately from testcross production.
- Use BPS trial/environment controls where possible.
- Custom allocation/mixed-model logic is likely needed for GCA, SCA, GxE, and environment correlations.
- Store MET hybrid results and `Training_RLines_GCA`.

`select_training_R_parents_by_GCA_index`

- Compute desired-gain index from GCA values across all traits.
- Select best 5 R lines as `R_Parents_Cycle_0`.
- Retain SCA as reportable/diagnostic.

`train_fixed_GS_model_from_GCA_index`

- Train one fixed GP model using genotypes and GCA/index response.
- Store as `GS_Model`.
- Do not retrain during the five recurrent cycles.
- Later implementation must decide whether to train directly on index values or train multi-trait GCA models and derive predicted index.

`make_recurrent_R_crosses`

- For each cycle, cross the current 5 parents with configurable 5-10 crosses.
- Store provenance by cycle.

`advance_recurrent_to_F5_and_genotype`

- Generate 50 total candidates across crosses.
- Advance by SSD/RGA to F5.
- Genotype all 50.
- Store `R_Candidates_F5_Cycle_i` at least until reporting and selection finish.

`predict_and_select_recurrent_R_parents`

- Cycles 1-4 only.
- Predict fixed-model index GEBV.
- Select 5.
- Store as `R_Parents_Cycle_i`.

`predict_and_select_final30_RLines`

- Cycle 5 only.
- Predict fixed-model index GEBV.
- Select 30 at F5.
- Store `R_Final30_F5`.

`advance_final30_F5_to_F6`

- Advance only the 30 selected lines to F6.
- Store `R_Final30_F6`.

`make_final_testcross_hybrids`

- Cross all 30 selected R lines to both fixed testers.
- Produce 60 candidate hybrid combinations.
- Store as `AYT_Hybrids`.

`run_partner_AYT`

- Run partner AYT on 60 hybrids for the same trait set.
- Store observed AYT data.
- Do not retrain the current model.

`release_hybrids_from_AYT`

- Release every hybrid whose observed AYT yield strictly outperforms the local commercial check.
- Store `Released_Hybrids`.
- Report other traits but do not gate release on them in this version.

## Deferred Implementation And Validation Work

Future implementation session should:

1. Convert this handoff into a BPS 0.2.0 scheme script.
2. Define all cfg fields explicitly.
3. Implement connected training testcross production and MET allocation at the v2 default scale of 800.
4. Verify GCA/SCA/GxE estimation approach.
5. Decide whether the GP target is direct index value or multi-trait GCA values.
6. Run `bp_check_cfg_requirements()`.
7. Create a minimal GP smoke test with `cfg$debug_GP` and `cfg$debug_GP_n`.
8. Validate event chronology, provenance, population sizes, costs, reporting, prediction, and release logic.

## Source References Consulted

The design used the v2 image plus selective references only:

- `Felicien_Variety_development_v2.jpeg`
- `.agents/skills/breeding-scheme-drafter/references/script-style-guide.md`
- `.agents/skills/breeding-scheme-drafter/references/discovery-guide.md`
- `BreedingProgramSimulator/inst/templates/bps-0.2.0/TEMPLATE_INDEX.md`
- `BreedingProgramSimulator/inst/templates/bps-0.2.0/scheme_phenotypic_multi_trait.R`
- `BreedingProgramSimulator/inst/templates/bps-0.2.0/scheme_two_part_gp_single_trait.R`
- `BreedingProgramSimulator/inst/templates/bps-0.2.0/scheme_wfrgs_single_trait.R`
- Narrow `rg` checks in `BreedingProgramSimulator/R` for BPS operations related to trials, genotyping, GP, synthetic traits, phenotype-allocation support, and progeny/testcross helpers.
