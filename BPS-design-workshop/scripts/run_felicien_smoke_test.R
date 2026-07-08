# Minimal smoke test for the Felicien hybrid rice scheme.
#
# Run from the workspace root:
#   source("scripts/run_felicien_smoke_test.R")
#
# This is intentionally tiny. It is for debugging event flow, cfg fields,
# source stages, GP integration, and release logic, not for biological results.

library(AlphaSimR)
library(BreedingProgramSimulator)
library(dplyr)

source("cfg/felicien_variety_development_cfg.R")
source("scripts/felicien_variety_development_scheme.R")

set.seed(101)

cfg_smoke <- cfg
cfg_smoke$n_training_crosses <- 5L
cfg_smoke$n_training_rlines <- 30L
cfg_smoke$training_met_n_locs <- 2L
cfg_smoke$training_met_max_plots_per_location <- 40L
cfg_smoke$training_met_env_means <- c(-0.2, 0.2)

cfg_smoke$n_recurrent_crosses <- 5L
cfg_smoke$n_recurrent_f5_candidates <- 12L
cfg_smoke$n_recurrent_parents <- 5L
cfg_smoke$recurrent_cross_duration_years <- 6 / 12
cfg_smoke$recurrent_ssd_duration_years <- 12 / 12
cfg_smoke$genotyping_duration_years <- 6 / 12
cfg_smoke$n_final_r_lines <- 6L

cfg_smoke$ayt_n_locs <- 2L
cfg_smoke$ayt_env_means <- c(-0.2, 0.2)
cfg_smoke$ayt_local_check_yield <- -Inf

cfg_smoke$debug_GP <- TRUE
cfg_smoke$debug_GP_n <- 20L

founder_pop <- AlphaSimR::quickHaplo(
  nInd = 80,
  nChr = 2,
  segSites = 60,
  inbred = TRUE
)

SP <- AlphaSimR::SimParam$new(founder_pop)
SP$restrSegSites(10, 30)
SP$addSnpChip(nSnpPerChr = 30)
SP$addTraitA(
  nQtlPerChr = 10,
  mean = rep(0, cfg_smoke$nTraits),
  var = cfg_smoke$traitVarG,
  name = cfg_smoke$trait_names
)
SP$setTrackPed(TRUE)

base_pop <- AlphaSimR::newPop(founder_pop, simParam = SP)
base_pop <- AlphaSimR::setPheno(
  pop = base_pop,
  varE = cfg_smoke$varE,
  reps = 1,
  traits = seq_len(cfg_smoke$nTraits),
  simParam = SP
)

r_source <- AlphaSimR::selectInd(
  pop = base_pop,
  nInd = 20,
  trait = AlphaSimR::selIndex,
  use = "pheno",
  b = cfg_smoke$initial_r_selection_weights,
  simParam = SP
)

tester_source <- AlphaSimR::selectInd(
  pop = base_pop,
  nInd = 10,
  trait = AlphaSimR::selIndex,
  use = "pheno",
  b = cfg_smoke$initial_tester_selection_weights,
  simParam = SP
)

state <- bp_init_state(
  SP = SP,
  dt = 1 / cfg_smoke$ticks_per_year,
  start_time = 0,
  sim = list(default_chip = cfg_smoke$snpChip)
)
state <- bp_register_synthetic_traits(state, cfg_smoke$synthetic_trait)

state <- put_stage_pop(
  state = state,
  pop = r_source,
  stage = "Parents",
  source_ids = "UNKNOWN",
  ready_in_years = 0,
  selection_strategy = "Smoke-test R-line source population",
  inherit_genotypes = TRUE,
  stream = "main"
)

state <- put_stage_pop(
  state = state,
  pop = tester_source,
  stage = "Parents",
  source_ids = "UNKNOWN",
  ready_in_years = 0,
  selection_strategy = "Smoke-test tester source population",
  inherit_genotypes = TRUE,
  stream = "tester"
)

smoke <- run_simulation(
  state = state,
  results_base = data.frame(Scenario = "felicien_smoke"),
  cfg = cfg_smoke,
  seed = 101,
  start_year = 0
)

cat("\nSmoke test completed.\n")
print(tail(smoke$state$cohorts[, c("cohort_id", "stage", "stream", "cycle_id", "n_ind", "available_tick")], 12))

recurrent_cycle_years <-
  cfg_smoke$recurrent_cross_duration_years +
  cfg_smoke$recurrent_ssd_duration_years +
  cfg_smoke$genotyping_duration_years
recurrent_cycle_months <- recurrent_cycle_years * 12

cat("\nConfigured recurrent GS cycle time for genetic gain:\n")
cat(sprintf(
  "- crossing: %.0f months\n",
  cfg_smoke$recurrent_cross_duration_years * 12
))
cat(sprintf(
  "- RGA/SSD to F5: %.0f months\n",
  cfg_smoke$recurrent_ssd_duration_years * 12
))
cat(sprintf(
  "- genotyping turnaround: %.0f months\n",
  cfg_smoke$genotyping_duration_years * 12
))
cat(sprintf(
  "- total L: %.0f months = %.2f years per recurrent GS cycle\n",
  recurrent_cycle_months,
  recurrent_cycle_years
))

parent_rows <- smoke$state$cohorts[
  smoke$state$cohorts$stage == "R_Parents" &
    smoke$state$cohorts$stream == "main",
  c("cohort_id", "cycle_id", "available_tick", "n_ind"),
  drop = FALSE
]
parent_rows$available_year <- parent_rows$available_tick * smoke$state$time$dt
parent_rows$delta_years_from_previous_parent_selection <- c(
  NA_real_,
  diff(parent_rows$available_year)
)
cat("\nObserved selected-parent timing by cycle:\n")
print(parent_rows)

gain_table <- parent_genetic_gain_table(smoke$state, cfg_smoke)
cat("\nRealized genetic gain in selected R parents:\n")
print(gain_table)

cat("\nIndex genetic gain summary:\n")
print(gain_table[, c(
  "cycle_id",
  "mean_index_gv",
  "delta_index_gv_per_cycle",
  "delta_index_gv_per_year"
), drop = FALSE])

if (nrow(gain_table) > 1L) {
  cumulative_gain <- tail(gain_table$mean_index_gv, 1L) - gain_table$mean_index_gv[[1L]]
  n_completed_cycles <- nrow(gain_table) - 1L
  mean_gain_per_cycle <- cumulative_gain / n_completed_cycles
  cat("\nCumulative recurrent-parent index gain:\n")
  cat(sprintf(
    "- cycle_0 to %s: %.4f index-GV units across %d recurrent cycles\n",
    tail(gain_table$cycle_id, 1L),
    cumulative_gain,
    n_completed_cycles
  ))
  cat(sprintf(
    "- mean gain per cycle: %.4f index-GV units/cycle\n",
    mean_gain_per_cycle
  ))
  cat(sprintf(
    "- mean gain per year at L = %.2f years: %.4f index-GV units/year\n",
    recurrent_cycle_years,
    mean_gain_per_cycle / recurrent_cycle_years
  ))
}

cat("\nEvent timeline:\n")
bp_print_event_timeline(smoke$state, collapse_year_patterns = FALSE)
