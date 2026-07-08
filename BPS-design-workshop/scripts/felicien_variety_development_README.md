# Felicien Variety Development BPS Draft

Files:

- `scripts/felicien_variety_development_scheme.R`
- `cfg/felicien_variety_development_cfg.R`
- `design_notes/felicien_variety_development_handoff.md`

This draft converts the approved design into an actual BPS scheme script and explicit cfg starter fields. It has not been run because R is not installed in this session.

## Implemented In The Scheme

- 5 fixed elite R-line founders.
- 2 fixed female testers.
- Training R-line development from 5-10 crosses to 800-1000 SSD/RGA lines.
- Sparse testcross generation with anchors for tester connectivity.
- Sparse partner MET allocation with a maximum plot count per location.
- GCA-style summaries from testcross hybrid phenotypes back to R lines.
- Desired-gain GCA index stored as a BPS synthetic phenotype.
- Fixed GP model trained from the training-cycle GCA index.
- Five recurrent GS cycles.
- Cycles 1-4 select 5 F5 parents by fixed-model index GEBV.
- Cycle 5 selects 30 F5 R lines, advances to F7, and crosses each to both testers.
- Final AYT and yield-check release rule.

## Needs Runtime Validation

- For cfg scanning, scan the scheme file and pass the cfg file only as
  `cfg_file`. Do not include the cfg file in `files`, or every intentional cfg
  assignment will be reported as "overwritten".
- Source the scheme without running it.
- Run `bp_check_cfg_requirements()` against the scheme and starter cfg.
- Verify the sparse MET helper against BPS phenotype logging expectations.
- The sparse MET helper currently stores aggregate entry phenotypes on the trial
  population and logs the allocation event; per-environment phenotype-log rows
  should be reviewed during validation.
- Verify the custom GCA/SCA approximation and decide whether to replace it with a mixed-model analysis.
- Verify fixed-model synthetic-trait GP training and prediction.
- Confirm that the setup state provides separate R-line and tester source stages.
- Run a minimal GP smoke test with very small population sizes before production runs.

Suggested cfg scan:

```r
bp_check_cfg_requirements(
  files = "scripts/felicien_variety_development_scheme.R",
  cfg_file = "cfg/felicien_variety_development_cfg.R",
  rewrite_file = FALSE
)
```

Suggested smoke test:

```r
source("scripts/run_felicien_smoke_test.R")
```

Suggested replicated comparisons:

```r
source("scripts/run_felicien_experiments.R")
```

By default this runs 20 replicates for each comparison scenario and writes:

- `results/felicien_experiments/scenario_design.csv`
- `results/felicien_experiments/replicate_summary.csv`
- `results/felicien_experiments/parent_gain_by_cycle.csv`
- `results/felicien_experiments/scenario_summary.csv`

For a quick runtime check from a terminal, run fewer replicates first:

```r
Rscript scripts/run_felicien_experiments.R 5
```
