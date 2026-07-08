# Starter cfg fields for scripts/felicien_variety_development_scheme.R.
#
# This file is intentionally explicit. Values here are placeholders for
# implementation and smoke-test setup; production values should be reviewed
# with the breeder before running experiments.

felicien_index <- function(x, weights) {
  as.numeric(as.matrix(x) %*% as.numeric(weights))
}

cfg <- list()

# Trait architecture -------------------------------------------------------
cfg$trait_names <- c(
  "yield",
  "blast",
  "hoja_blanca",
  "Bulkhoderia",
  "milling_yield",
  "milling_quality",
  "white_center"
)
cfg$nTraits <- length(cfg$trait_names)
cfg$traitNames <- cfg$trait_names

# Positive weights mean higher is preferred. Use negative weights for traits
# where lower observed values are preferred.
cfg$gca_index_weights <- c(
  yield = 1.00,
  blast = -0.30,
  hoja_blanca = -0.30,
  Bulkhoderia = -0.30,
  milling_yield = 0.35,
  milling_quality = 0.35,
  white_center = -0.20
)
cfg$synthetic_trait_name <- "felicien_gca_index"
cfg$synthetic_trait <- bp_synthetic_trait(
  name = cfg$synthetic_trait_name,
  traits = seq_len(cfg$nTraits),
  fun = felicien_index,
  args = list(weights = cfg$gca_index_weights),
  linear = TRUE,
  missing_component = "zero_if_all_missing"
)

# Initial source stages ----------------------------------------------------
# The implementation expects the setup state to contain an R-line source stage
# and a tester source stage. These names can be changed to match Create_sim.
cfg$initial_r_source_stage <- "Parents"
cfg$initial_r_source_stream <- "main"
cfg$initial_tester_source_stage <- "Parents"
cfg$initial_tester_source_stream <- "tester"
cfg$n_r_founders <- 5L
cfg$n_testers <- 2L
cfg$initial_r_selection_use <- "pheno"
cfg$initial_tester_selection_use <- "pheno"
cfg$initial_r_selection_weights <- cfg$gca_index_weights
cfg$initial_tester_selection_weights <- cfg$gca_index_weights

# Genotyping and GP --------------------------------------------------------
cfg$snpChip <- 1L
cfg$genotyping_duration_years <- 0.5
cfg$genotype_cost_per_sample <- 15
cfg$gs_model_id <- "felicien_fixed_gca_index_model"
cfg$training_model_lookback_years <- 20
cfg$debug_GP <- TRUE
cfg$debug_GP_n <- 100L

# Training cycle -----------------------------------------------------------
cfg$n_training_crosses <- 10L
cfg$n_training_rlines <- 800L
cfg$training_cross_duration_years <- 6 / 12
cfg$training_ssd_selfing_generations <- 4L
cfg$training_ssd_duration_years <- 12 / 12
cfg$training_testcross_duration_years <- 5 / 12
cfg$training_testcross_anchor_fraction <- 0.10

cfg$training_cross_cost_per_cross <- 25
cfg$training_line_development_cost_per_line <- 1
cfg$training_seed_increase_cost_per_line <- 1
cfg$training_testcross_cost_per_cross <- 20

# Partner MET: two sparse layers are represented in the scheme:
# 1. sparse testcrossing across two testers
# 2. sparse phenotyping allocation across locations
cfg$training_met_n_locs <- 5L
cfg$training_met_reps <- 1L
cfg$training_met_max_plots_per_location <- 200L
cfg$training_met_duration_years <- 1.0
cfg$training_met_cost_per_plot <- 30
cfg$training_met_env_means <- c(-0.5, -0.2, 0.0, 0.2, 0.5)
cfg$training_met_env_mean_sd <- 0.25
cfg$training_met_env_year_sd <- 0.20

# Recurrent GS cycles ------------------------------------------------------
cfg$n_recurrent_cycles <- 5L
cfg$n_recurrent_crosses <- 10L
cfg$n_recurrent_f5_candidates <- 50L
cfg$n_recurrent_parents <- 5L
cfg$recurrent_cross_duration_years <- 6 / 12
cfg$recurrent_ssd_selfing_generations <- 4L
cfg$recurrent_ssd_duration_years <- 12 / 12

cfg$recurrent_cross_cost_per_cross <- 25
cfg$recurrent_line_development_cost_per_line <- 1
cfg$recurrent_seed_increase_cost_per_line <- 1

# Final product pipeline ---------------------------------------------------
cfg$n_final_r_lines <- 30L
cfg$final_f5_to_f7_selfing_generations <- 2L
cfg$final_f5_to_f7_duration_years <- 6 / 12
cfg$final_line_development_cost_per_line <- 1
cfg$final_testcross_duration_years <- 5 / 12
cfg$final_testcross_cost_per_cross <- 20

cfg$ayt_n_locs <- 5L
cfg$ayt_reps <- 1L
cfg$ayt_duration_years <- 1.0
cfg$ayt_cost_per_plot <- 35
cfg$ayt_env_means <- c(-0.5, -0.2, 0.0, 0.2, 0.5)
cfg$ayt_env_mean_sd <- 0.25
cfg$ayt_env_year_sd <- 0.20
cfg$ayt_local_check_yield <- 0

# Genetic/residual placeholders ------------------------------------------
# Create_sim should replace these with calibrated values.
cfg$traitVarG <- rep(1, cfg$nTraits)
cfg$traitH2 <- rep(0.35, cfg$nTraits)
cfg$traitCorE <- diag(cfg$nTraits)
cfg$varE <- bp_make_varE(
  h2 = cfg$traitH2,
  varG = cfg$traitVarG,
  corE = cfg$traitCorE,
  trait_names = cfg$trait_names
)

# Reporting / orchestration -----------------------------------------------
cfg$ticks_per_year <- 12L
cfg$nYearsFuture <- 1L
