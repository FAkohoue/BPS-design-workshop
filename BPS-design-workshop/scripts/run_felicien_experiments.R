# Production-style replicated comparisons for the Felicien hybrid rice scheme.
#
# Run from the workspace root in RStudio:
#   source("scripts/run_felicien_experiments.R")
#
# Or from a terminal:
#   Rscript scripts/run_felicien_experiments.R 20
#
# The optional first command-line argument is the number of replicates per
# scenario. Start with 5-10, then move toward 50-100 once runtime is known.

library(AlphaSimR)
library(BreedingProgramSimulator)
library(dplyr)

source("cfg/felicien_variety_development_cfg.R")
source("scripts/felicien_variety_development_scheme.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1L && nzchar(args[[1]])) as.integer(args[[1]]) else 100L
if (is.na(n_reps) || n_reps < 1L) stop("n_reps must be a positive integer.", call. = FALSE)

results_dir <- file.path("results", "felicien_experiments")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

scenario_design <- data.frame(
  ScenarioID = 1:5,
  Scenario = c(
    "baseline_800_met12_ayt12",
    "fast_partner_data_800_met6_ayt6",
    "larger_training_900_met12_ayt12",
    "fewer_crosses_5_train800_met12_ayt12",
    "higher_capacity_1000_cap250_met12_ayt12"
  ),
  n_training_crosses = c(10L, 10L, 10L, 5L, 10L),
  n_training_rlines = c(800L, 800L, 900L, 800L, 1000L),
  training_met_max_plots_per_location = c(200L, 200L, 200L, 200L, 250L),
  training_met_duration_years = c(12 / 12, 6 / 12, 12 / 12, 12 / 12, 12 / 12),
  ayt_duration_years = c(12 / 12, 6 / 12, 12 / 12, 12 / 12, 12 / 12),
  recurrent_cross_duration_years = c(6 / 12, 6 / 12, 6 / 12, 6 / 12, 6 / 12),
  recurrent_ssd_duration_years = c(12 / 12, 12 / 12, 12 / 12, 12 / 12, 12 / 12),
  genotyping_duration_years = c(6 / 12, 6 / 12, 6 / 12, 6 / 12, 6 / 12),
  stringsAsFactors = FALSE
)

make_experiment_cfg <- function(base_cfg, scenario_row) {
  cfg_run <- base_cfg
  for (nm in names(scenario_row)) {
    if (nm %in% c("ScenarioID", "Scenario")) next
    cfg_run[[nm]] <- scenario_row[[nm]]
  }
  cfg_run$ScenarioID <- scenario_row$ScenarioID
  cfg_run$Scenario <- scenario_row$Scenario
  cfg_run$debug_GP <- FALSE
  cfg_run
}

make_initial_state <- function(cfg_run, seed) {
  set.seed(seed)

  founder_pop <- AlphaSimR::quickHaplo(
    nInd = 300,
    nChr = 12,
    segSites = 600,
    inbred = TRUE
  )

  SP <- AlphaSimR::SimParam$new(founder_pop)
  SP$restrSegSites(100, 500)
  SP$addSnpChip(nSnpPerChr = 500)
  SP$addTraitA(
    nQtlPerChr = 100,
    mean = rep(0, cfg_run$nTraits),
    var = cfg_run$traitVarG,
    name = cfg_run$trait_names
  )
  SP$setTrackPed(TRUE)

  base_pop <- AlphaSimR::newPop(founder_pop, simParam = SP)
  base_pop <- AlphaSimR::setPheno(
    pop = base_pop,
    varE = cfg_run$varE,
    reps = 1,
    traits = seq_len(cfg_run$nTraits),
    simParam = SP
  )

  r_source <- AlphaSimR::selectInd(
    pop = base_pop,
    nInd = 80,
    trait = AlphaSimR::selIndex,
    use = "pheno",
    b = cfg_run$initial_r_selection_weights,
    simParam = SP
  )

  tester_source <- AlphaSimR::selectInd(
    pop = base_pop,
    nInd = 20,
    trait = AlphaSimR::selIndex,
    use = "pheno",
    b = cfg_run$initial_tester_selection_weights,
    simParam = SP
  )

  state <- bp_init_state(
    SP = SP,
    dt = 1 / cfg_run$ticks_per_year,
    start_time = 0,
    sim = list(default_chip = cfg_run$snpChip)
  )
  state <- bp_register_synthetic_traits(state, cfg_run$synthetic_trait)

  state <- put_stage_pop(
    state = state,
    pop = r_source,
    stage = "Parents",
    source_ids = "UNKNOWN",
    ready_in_years = 0,
    selection_strategy = "Experiment R-line source population",
    inherit_genotypes = TRUE,
    stream = "main"
  )

  put_stage_pop(
    state = state,
    pop = tester_source,
    stage = "Parents",
    source_ids = "UNKNOWN",
    ready_in_years = 0,
    selection_strategy = "Experiment tester source population",
    inherit_genotypes = TRUE,
    stream = "tester"
  )
}

summarize_replicate <- function(sim, cfg_run, scenario_row, replicate_id, seed) {
  gain <- parent_genetic_gain_table(sim$state, cfg_run)
  latest_gain <- latest_parent_gain_metrics(sim$state, cfg_run)
  released <- select_latest_available(
    state = sim$state,
    stage = "Released_Hybrids",
    stream = "main",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )
  ayt <- select_latest_available(
    state = sim$state,
    stage = "AYT_Results",
    stream = "main",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )

  data.frame(
    ScenarioID = scenario_row$ScenarioID,
    Scenario = scenario_row$Scenario,
    ReplicateID = replicate_id,
    Seed = seed,
    recurrent_cycle_time_years = recurrent_cycle_time_years(cfg_run),
    n_training_rlines = cfg_run$n_training_rlines,
    training_met_max_plots_per_location = cfg_run$training_met_max_plots_per_location,
    n_recurrent_f5_candidates = cfg_run$n_recurrent_f5_candidates,
    n_final_r_lines = cfg_run$n_final_r_lines,
    training_met_duration_years = cfg_run$training_met_duration_years,
    ayt_duration_years = cfg_run$ayt_duration_years,
    n_parent_cycles = nrow(gain),
    latest_parent_mean_index_gv = latest_gain[["latest_parent_mean_index_gv"]],
    cumulative_parent_index_gain = latest_gain[["cumulative_parent_index_gain"]],
    mean_parent_index_gain_per_cycle = latest_gain[["mean_parent_index_gain_per_cycle"]],
    mean_parent_index_gain_per_year = latest_gain[["mean_parent_index_gain_per_year"]],
    n_released_hybrids = if (is.null(released)) 0L else nInd(released$pop),
    mean_released_yield_gv = if (is.null(released)) NA_real_ else mean(released$pop@gv[, "yield"], na.rm = TRUE),
    mean_ayt_yield_pheno = if (is.null(ayt)) NA_real_ else mean(ayt$pop@pheno[, match("yield", cfg_run$trait_names)], na.rm = TRUE),
    final_year = sim$state$time$t,
    stringsAsFactors = FALSE
  )
}

se_value <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

all_rep_summaries <- list()
all_parent_gain <- list()

for (s in seq_len(nrow(scenario_design))) {
  scenario_row <- scenario_design[s, , drop = FALSE]
  cfg_run <- make_experiment_cfg(cfg, scenario_row)

  cat(sprintf("\nScenario %d/%d: %s\n", s, nrow(scenario_design), scenario_row$Scenario))

  for (replicate_id in seq_len(n_reps)) {
    seed <- scenario_row$ScenarioID * 10000L + replicate_id
    cat(sprintf("  replicate %d/%d seed=%d\n", replicate_id, n_reps, seed))

    state <- make_initial_state(cfg_run, seed = seed)
    sim <- run_simulation(
      state = state,
      results_base = data.frame(
        ScenarioID = scenario_row$ScenarioID,
        Scenario = scenario_row$Scenario,
        ReplicateID = replicate_id,
        Seed = seed
      ),
      cfg = cfg_run,
      seed = seed,
      start_year = 0
    )

    rep_summary <- summarize_replicate(
      sim = sim,
      cfg_run = cfg_run,
      scenario_row = scenario_row,
      replicate_id = replicate_id,
      seed = seed
    )
    gain <- parent_genetic_gain_table(sim$state, cfg_run)
    if (nrow(gain)) {
      gain$ScenarioID <- scenario_row$ScenarioID
      gain$Scenario <- scenario_row$Scenario
      gain$ReplicateID <- replicate_id
      gain$Seed <- seed
    }

    all_rep_summaries[[length(all_rep_summaries) + 1L]] <- rep_summary
    all_parent_gain[[length(all_parent_gain) + 1L]] <- gain
  }
}

replicate_summary <- bind_rows(all_rep_summaries)
parent_gain <- bind_rows(all_parent_gain)

scenario_summary <- replicate_summary %>%
  group_by(ScenarioID, Scenario) %>%
  summarise(
    n_reps = n(),
    n_valid_gain_reps = sum(is.finite(mean_parent_index_gain_per_cycle)),
    recurrent_cycle_time_years = mean(recurrent_cycle_time_years),
    se_parent_index_gain_per_cycle = se_value(mean_parent_index_gain_per_cycle),
    mean_parent_index_gain_per_cycle = mean(mean_parent_index_gain_per_cycle, na.rm = TRUE),
    se_parent_index_gain_per_year = se_value(mean_parent_index_gain_per_year),
    mean_parent_index_gain_per_year = mean(mean_parent_index_gain_per_year, na.rm = TRUE),
    mean_cumulative_parent_index_gain = mean(cumulative_parent_index_gain, na.rm = TRUE),
    mean_n_released_hybrids = mean(n_released_hybrids, na.rm = TRUE),
    mean_final_year = mean(final_year, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  scenario_design,
  file = file.path(results_dir, "scenario_design.csv"),
  row.names = FALSE
)
write.csv(
  replicate_summary,
  file = file.path(results_dir, "replicate_summary.csv"),
  row.names = FALSE
)
write.csv(
  parent_gain,
  file = file.path(results_dir, "parent_gain_by_cycle.csv"),
  row.names = FALSE
)
write.csv(
  scenario_summary,
  file = file.path(results_dir, "scenario_summary.csv"),
  row.names = FALSE
)

cat("\nExperiment summaries written to:\n")
cat(sprintf("- %s\n", normalizePath(results_dir, winslash = "/", mustWork = FALSE)))
cat("\nScenario summary:\n")
print(scenario_summary)

cat("\nCompact gain comparison with standard errors:\n")
print(as.data.frame(scenario_summary[, c(
  "ScenarioID",
  "Scenario",
  "n_valid_gain_reps",
  "mean_parent_index_gain_per_cycle",
  "se_parent_index_gain_per_cycle",
  "mean_parent_index_gain_per_year",
  "se_parent_index_gain_per_year",
  "mean_cumulative_parent_index_gain",
  "mean_final_year"
)]), row.names = FALSE)
