# Setup -------------------------------------------------------------------
# Felicien hybrid rice variety-development scheme for BPS 0.2.0.
#
# Biological flow:
#   5 related elite R founders + 2 fixed female testers
#   -> training R-line family
#   -> sparse testcross MET
#   -> GCA desired-gain index
#   -> fixed GS model
#   -> 5 recurrent GS cycles
#   -> 60 final hybrids
#   -> partner AYT
#   -> release hybrids exceeding local yield check
#
# This is an implementation draft from
# design_notes/felicien_variety_development_handoff.md.
# R runtime validation is still required.

library(AlphaSimR)
library(BreedingProgramSimulator)
library(dplyr)

# Helper Utilities ---------------------------------------------------------

cost_value <- function(state, event) {
  if (is.null(state$cost_log) || !nrow(state$cost_log)) return(0)
  vals <- state$cost_log$total_cost[state$cost_log$event == event]
  if (!length(vals)) return(0)
  sum(vals, na.rm = TRUE)
}

total_cost_value <- function(state) {
  if (is.null(state$cost_log) || !nrow(state$cost_log)) return(0)
  sum(state$cost_log$total_cost, na.rm = TRUE)
}

trait_indices <- function(cfg) {
  seq_len(as.integer(cfg$nTraits))
}

make_pair_cross_plan <- function(n_parents, n_crosses, allow_selfing = FALSE) {
  n_parents <- as.integer(n_parents)
  n_crosses <- as.integer(n_crosses)
  if (n_parents < 2L) stop("At least two parents are required for crossing.", call. = FALSE)
  if (n_crosses < 1L) stop("n_crosses must be positive.", call. = FALSE)

  pair_pool <- if (isTRUE(allow_selfing)) {
    as.matrix(expand.grid(mother = seq_len(n_parents), father = seq_len(n_parents)))
  } else {
    t(utils::combn(seq_len(n_parents), 2L))
  }
  pair_pool[sample(seq_len(nrow(pair_pool)), n_crosses, replace = TRUE), , drop = FALSE]
}

make_crosses_from_parents <- function(parents, n_crosses, n_progeny_per_cross, simParam) {
  cross_plan <- make_pair_cross_plan(
    n_parents = nInd(parents),
    n_crosses = n_crosses,
    allow_selfing = FALSE
  )
  AlphaSimR::makeCross(
    pop = parents,
    crossPlan = cross_plan,
    nProgeny = as.integer(n_progeny_per_cross),
    simParam = simParam
  )
}

advance_by_ssd <- function(f1_pop, n_lines, n_selfing_generations, simParam) {
  n_lines <- as.integer(n_lines)
  n_selfing_generations <- as.integer(n_selfing_generations)
  if (n_lines < 1L) stop("n_lines must be positive.", call. = FALSE)
  if (n_selfing_generations < 1L) {
    stop("n_selfing_generations must be at least one.", call. = FALSE)
  }

  sampled_f1 <- sample(seq_len(nInd(f1_pop)), n_lines, replace = TRUE)
  out <- AlphaSimR::self(
    pop = f1_pop[sampled_f1],
    nProgeny = 1L,
    simParam = simParam
  )
  if (n_selfing_generations > 1L) {
    for (gen in seq_len(n_selfing_generations)[-1L]) {
      out <- AlphaSimR::self(pop = out, nProgeny = 1L, simParam = simParam)
    }
  }
  out
}

make_sparse_testcross_plan <- function(n_rlines, n_testers, anchor_fraction) {
  n_rlines <- as.integer(n_rlines)
  n_testers <- as.integer(n_testers)
  if (n_rlines < 1L || n_testers < 1L) stop("Invalid testcross dimensions.", call. = FALSE)

  base <- data.frame(
    r_index = seq_len(n_rlines),
    tester_index = ((seq_len(n_rlines) - 1L) %% n_testers) + 1L
  )
  n_anchor <- min(n_rlines, max(n_testers, ceiling(n_rlines * anchor_fraction)))
  anchor_lines <- seq_len(n_anchor)
  anchors <- expand.grid(
    r_index = anchor_lines,
    tester_index = seq_len(n_testers)
  )
  unique(rbind(base, anchors))
}

make_full_testcross_plan <- function(n_rlines, n_testers) {
  expand.grid(
    r_index = seq_len(as.integer(n_rlines)),
    tester_index = seq_len(as.integer(n_testers))
  )
}

make_testcross_pop <- function(rlines, testers, testcross_plan, simParam) {
  combined <- c(rlines, testers)
  n_r <- nInd(rlines)
  cross_plan <- cbind(
    as.integer(testcross_plan$r_index),
    n_r + as.integer(testcross_plan$tester_index)
  )
  hybrids <- AlphaSimR::makeCross(
    pop = combined,
    crossPlan = cross_plan,
    nProgeny = 1L,
    simParam = simParam
  )
  hybrids@misc[["r_line_id"]] <- rlines@id[as.integer(testcross_plan$r_index)]
  hybrids@misc[["tester_id"]] <- testers@id[as.integer(testcross_plan$tester_index)]
  hybrids
}

make_sparse_environment_allocation <- function(n_entries, n_locations, max_plots_per_location) {
  n_entries <- as.integer(n_entries)
  n_locations <- as.integer(n_locations)
  max_plots_per_location <- as.integer(max_plots_per_location)
  if (n_entries < 1L || n_locations < 1L) stop("Invalid sparse trial dimensions.", call. = FALSE)
  if (max_plots_per_location < 1L) stop("max_plots_per_location must be positive.", call. = FALSE)
  if (n_entries > n_locations * max_plots_per_location) {
    stop("Sparse MET capacity is smaller than the number of entries.", call. = FALSE)
  }

  allocation <- vector("list", n_locations)
  for (loc in seq_len(n_locations)) allocation[[loc]] <- integer(0)

  order_entries <- sample(seq_len(n_entries))
  for (i in seq_along(order_entries)) {
    loc <- ((i - 1L) %% n_locations) + 1L
    allocation[[loc]] <- c(allocation[[loc]], order_entries[[i]])
  }

  anchors <- seq_len(min(n_entries, n_locations))
  for (loc in seq_len(n_locations)) {
    free_slots <- max_plots_per_location - length(allocation[[loc]])
    if (free_slots > 0L) {
      add_anchors <- setdiff(anchors, allocation[[loc]])
      allocation[[loc]] <- unique(c(allocation[[loc]], head(add_anchors, free_slots)))
    }
  }
  allocation
}

resolve_trial_env_local <- function(env_means, env_mean_sd, env_year_sd, n_loc) {
  env_means <- as.numeric(env_means)
  n_loc <- as.integer(n_loc)
  if (length(env_means) != n_loc) {
    stop("Environment means length must equal number of locations.", call. = FALSE)
  }
  loc_dev <- stats::rnorm(n_loc, mean = 0, sd = as.numeric(env_mean_sd))
  year_eff <- stats::rnorm(1L, mean = 0, sd = as.numeric(env_year_sd))
  env_means + loc_dev + year_eff
}

env_p_from_latent_local <- function(latent_env, n_traits) {
  rep(stats::pnorm(as.numeric(latent_env)), as.integer(n_traits))
}

last_cohort_id_local <- function(state) {
  tail(as.character(state$cohorts$cohort_id), 1L)
}

assign_trial_pheno <- function(pop, traits, pheno_matrix) {
  full <- matrix(NA_real_, nrow = nInd(pop), ncol = length(traits))
  full[, seq_along(traits)] <- pheno_matrix
  colnames(full) <- colnames(pop@gv)[as.integer(traits)]
  pop@pheno <- full
  pop
}

run_sparse_partner_trial <- function(
  state,
  pop,
  output_stage,
  input_cohorts,
  cfg,
  trial_cfg_prefix,
  selection_strategy,
  stream = "main",
  cycle_id = "cycle_1"
) {
  traits <- trait_indices(cfg)
  n_loc <- as.integer(cfg[[paste0(trial_cfg_prefix, "_n_locs")]])
  reps <- as.integer(cfg[[paste0(trial_cfg_prefix, "_reps")]])
  max_plots <- as.integer(cfg[[paste0(trial_cfg_prefix, "_max_plots_per_location")]])
  duration_years <- as.numeric(cfg[[paste0(trial_cfg_prefix, "_duration_years")]])
  cost_per_plot <- as.numeric(cfg[[paste0(trial_cfg_prefix, "_cost_per_plot")]])
  env_means <- cfg[[paste0(trial_cfg_prefix, "_env_means")]]
  env_mean_sd <- as.numeric(cfg[[paste0(trial_cfg_prefix, "_env_mean_sd")]])
  env_year_sd <- as.numeric(cfg[[paste0(trial_cfg_prefix, "_env_year_sd")]])

  allocation <- make_sparse_environment_allocation(
    n_entries = nInd(pop),
    n_locations = n_loc,
    max_plots_per_location = max_plots
  )

  z_env <- resolve_trial_env_local(
    env_means = env_means,
    env_mean_sd = env_mean_sd,
    env_year_sd = env_year_sd,
    n_loc = n_loc
  )

  pheno_sum <- matrix(0, nrow = nInd(pop), ncol = length(traits))
  pheno_n <- integer(nInd(pop))

  for (loc in seq_len(n_loc)) {
    idx <- allocation[[loc]]
    p_env <- env_p_from_latent_local(latent_env = z_env[[loc]], n_traits = length(traits))
    ph_loc <- AlphaSimR::setPheno(
      pop = pop[idx],
      varE = cfg$varE[traits, traits],
      reps = reps,
      traits = traits,
      p = p_env,
      onlyPheno = TRUE,
      simParam = state$sim$SP
    )
    if (is.null(dim(ph_loc))) ph_loc <- matrix(ph_loc, ncol = length(traits))
    pheno_sum[idx, ] <- pheno_sum[idx, , drop = FALSE] + ph_loc
    pheno_n[idx] <- pheno_n[idx] + 1L
  }

  pheno_mean <- pheno_sum / pheno_n
  pop_trial <- assign_trial_pheno(pop, traits, pheno_mean)

  state <- put_stage_pop(
    state = state,
    pop = pop_trial,
    stage = output_stage,
    source = input_cohorts,
    ready_in_years = duration_years,
    stream = stream,
    cycle_id = cycle_id,
    inherit_genotypes = TRUE,
    selection_strategy = selection_strategy,
    cost_per_unit = cost_per_plot,
    cost_units = sum(lengths(allocation)) * reps,
    cost_event = "phenotype_trial",
    cost_unit = "plot"
  )

  bp_log_event(
    state = state,
    fn = "run_sparse_partner_trial",
    event_type = "sparse_phenotyping",
    stage = output_stage,
    source_ids = input_cohorts,
    output_id = last_cohort_id_local(state),
    event_string = sprintf(
      "Sparse partner trial for %s used %d locations, max %d plots/location, and %d total plot records before reps.",
      output_stage,
      n_loc,
      max_plots,
      sum(lengths(allocation))
    ),
    template_string = "Sparse partner MET/AYT allocation",
    details = list(
      allocation = allocation,
      n_loc = n_loc,
      max_plots_per_location = max_plots,
      reps = reps
    )
  )
}

attach_gca_sca_summaries <- function(rlines, hybrids, cfg) {
  r_ids <- as.integer(hybrids@misc[["r_line_id"]])
  tester_ids <- as.integer(hybrids@misc[["tester_id"]])
  ph <- as.matrix(hybrids@pheno)
  colnames(ph) <- cfg$trait_names

  gca <- matrix(NA_real_, nrow = nInd(rlines), ncol = ncol(ph))
  colnames(gca) <- paste0("gca_", cfg$trait_names)
  for (i in seq_len(nInd(rlines))) {
    idx <- r_ids == rlines@id[[i]]
    if (any(idx)) gca[i, ] <- colMeans(ph[idx, , drop = FALSE], na.rm = TRUE)
  }

  trait_center <- colMeans(gca, na.rm = TRUE)
  trait_scale <- apply(gca, 2L, stats::sd, na.rm = TRUE)
  trait_scale[!is.finite(trait_scale) | trait_scale == 0] <- 1
  gca_z <- sweep(sweep(gca, 2L, trait_center, "-"), 2L, trait_scale, "/")
  index <- as.numeric(gca_z %*% as.numeric(cfg$gca_index_weights))

  rlines@misc[["felicien_gca"]] <- gca
  rlines@misc[["felicien_gca_index"]] <- index
  rlines <- bp_set_synthetic_values(
    pop = rlines,
    name = cfg$synthetic_trait_name,
    values = index,
    type = "pheno"
  )

  sca <- ph
  for (j in seq_len(ncol(ph))) {
    r_mean <- tapply(ph[, j], r_ids, mean, na.rm = TRUE)
    t_mean <- tapply(ph[, j], tester_ids, mean, na.rm = TRUE)
    overall <- mean(ph[, j], na.rm = TRUE)
    sca[, j] <- ph[, j] - r_mean[as.character(r_ids)] - t_mean[as.character(tester_ids)] + overall
  }
  colnames(sca) <- paste0("sca_", cfg$trait_names)
  hybrids@misc[["felicien_sca"]] <- sca

  list(rlines = rlines, hybrids = hybrids)
}

score_latest_candidates_with_fixed_model <- function(state, src, cfg, stage_label) {
  state <- run_predict_ebv(
    state,
    list(
      cohort_ids = src$source_ids,
      model_id = cfg$gs_model_id,
      chip = cfg$snpChip,
      require_genotyped = TRUE,
      synthetic_trait = cfg$synthetic_trait_name
    )
  )
  scored <- select_latest_available(
    state = state,
    stage = stage_label,
    stream = "main",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )
  list(state = state, scored = scored)
}

recurrent_cycle_time_years <- function(cfg) {
  as.numeric(cfg$recurrent_cross_duration_years) +
    as.numeric(cfg$recurrent_ssd_duration_years) +
    as.numeric(cfg$genotyping_duration_years)
}

parent_genetic_gain_table <- function(state, cfg) {
  rows <- state$cohorts[
    state$cohorts$stage == "R_Parents" &
      state$cohorts$stream == "main",
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) return(data.frame())
  rows <- rows[order(rows$available_tick, rows$cohort_id), , drop = FALSE]

  out <- lapply(seq_len(nrow(rows)), function(i) {
    cid <- as.character(rows$cohort_id[[i]])
    pop <- state$pops[[cid]]
    index_gv <- bp_synthetic_values(pop, cfg$synthetic_trait, use = "gv")
    trait_gv <- colMeans(pop@gv[, seq_len(cfg$nTraits), drop = FALSE], na.rm = TRUE)
    names(trait_gv) <- paste0("mean_gv_", cfg$trait_names)

    data.frame(
      cohort_id = cid,
      cycle_id = as.character(rows$cycle_id[[i]]),
      available_year = as.numeric(rows$available_tick[[i]]) * state$time$dt,
      n_ind = as.integer(rows$n_ind[[i]]),
      mean_index_gv = mean(index_gv, na.rm = TRUE),
      as.list(trait_gv),
      check.names = FALSE
    )
  })
  out <- bind_rows(out)

  out$delta_index_gv_per_cycle <- c(NA_real_, diff(out$mean_index_gv))
  out$delta_index_gv_per_year <- out$delta_index_gv_per_cycle / recurrent_cycle_time_years(cfg)

  trait_cols <- paste0("mean_gv_", cfg$trait_names)
  for (trait_col in trait_cols) {
    delta_name <- sub("^mean_gv_", "delta_gv_", trait_col)
    out[[paste0(delta_name, "_per_cycle")]] <- c(NA_real_, diff(out[[trait_col]]))
    out[[paste0(delta_name, "_per_year")]] <-
      out[[paste0(delta_name, "_per_cycle")]] / recurrent_cycle_time_years(cfg)
  }

  out
}

latest_parent_gain_metrics <- function(state, cfg) {
  gain <- parent_genetic_gain_table(state, cfg)
  if (!nrow(gain)) {
    return(c(
      recurrent_cycle_time_years = recurrent_cycle_time_years(cfg),
      latest_parent_mean_index_gv = NA_real_,
      latest_parent_delta_index_gv_per_cycle = NA_real_,
      latest_parent_delta_index_gv_per_year = NA_real_,
      cumulative_parent_index_gain = NA_real_,
      mean_parent_index_gain_per_cycle = NA_real_,
      mean_parent_index_gain_per_year = NA_real_
    ))
  }
  latest <- gain[nrow(gain), , drop = FALSE]
  n_completed_cycles <- max(0L, nrow(gain) - 1L)
  cumulative_gain <- if (n_completed_cycles > 0L) {
    latest$mean_index_gv - gain$mean_index_gv[[1L]]
  } else {
    NA_real_
  }
  mean_gain_per_cycle <- if (n_completed_cycles > 0L) cumulative_gain / n_completed_cycles else NA_real_
  c(
    recurrent_cycle_time_years = recurrent_cycle_time_years(cfg),
    latest_parent_mean_index_gv = latest$mean_index_gv,
    latest_parent_delta_index_gv_per_cycle = latest$delta_index_gv_per_cycle,
    latest_parent_delta_index_gv_per_year = latest$delta_index_gv_per_year,
    cumulative_parent_index_gain = cumulative_gain,
    mean_parent_index_gain_per_cycle = mean_gain_per_cycle,
    mean_parent_index_gain_per_year = mean_gain_per_cycle / recurrent_cycle_time_years(cfg)
  )
}

# Reporting ----------------------------------------------------------------

record_yearly_outputs <- function(results, results_base, state, year, cfg) {
  all_traits <- trait_indices(cfg)
  results_year <- c(
    year = year,
    bp_report_stage_metrics(
      state = state,
      stage = "Training_RLines_GCA",
      stream = "main",
      metrics = c(meanTraining = "meanG", varTraining = "varG"),
      traits = all_traits,
      append_trait = "always"
    ),
    bp_report_stage_metrics(
      state = state,
      stage = "R_Candidates_F5",
      stream = "main",
      metrics = c(meanF5 = "meanG", varF5 = "varG", accF5 = "accEBV"),
      traits = all_traits,
      append_trait = "always"
    ),
    bp_report_stage_metrics(
      state = state,
      stage = "AYT_Hybrids",
      stream = "main",
      metrics = c(meanAYT = "meanG", H2_AYT = "H2"),
      traits = all_traits,
      append_trait = "always"
    ),
    bp_report_stage_metrics(
      state = state,
      stage = "Released_Hybrids",
      stream = "main",
      metrics = c(meanReleased = "meanG"),
      traits = all_traits,
      append_trait = "always"
    ),
    crossing_cost = cost_value(state, "crossing"),
    line_development_cost = cost_value(state, "line_development"),
    seed_increase_cost = cost_value(state, "seed_increase"),
    phenotype_cost = cost_value(state, "phenotype_trial"),
    genotype_cost = cost_value(state, "genotyping"),
    latest_parent_gain_metrics(state, cfg),
    total_cost = total_cost_value(state)
  )
  results_year[lengths(results_year) == 0] <- matrix(NA)
  bind_rows(results, data.frame(results_base, results_year))
}

# Event Verbs --------------------------------------------------------------

initialize_R_founders <- function(state, cfg, year) {
  existing <- select_latest_available(
    state = state,
    stage = "R_Founders",
    stream = "main",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )
  if (!is.null(existing)) return(state)

  source <- select_latest_available(
    state = state,
    stage = cfg$initial_r_source_stage,
    stream = cfg$initial_r_source_stream,
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )
  chk <- bp_skip_if_no_input(state, source, cfg, event_name = "initialize_R_founders")
  if (chk$skip) return(chk$state)

  founders <- AlphaSimR::selectInd(
    pop = source$pop,
    nInd = as.integer(cfg$n_r_founders),
    trait = AlphaSimR::selIndex,
    use = cfg$initial_r_selection_use,
    b = cfg$initial_r_selection_weights,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = founders,
    stage = "R_Founders",
    source = source,
    ready_in_years = 0,
    inherit_genotypes = TRUE,
    selection_strategy = "Fixed elite R founders selected from one outstanding family",
    stream = "main"
  )
}

initialize_fixed_testers <- function(state, cfg, year) {
  existing <- select_latest_available(
    state = state,
    stage = "Female_Testers",
    stream = "tester",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )
  if (!is.null(existing)) return(state)

  source <- select_latest_available(
    state = state,
    stage = cfg$initial_tester_source_stage,
    stream = cfg$initial_tester_source_stream,
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )
  chk <- bp_skip_if_no_input(state, source, cfg, event_name = "initialize_fixed_testers")
  if (chk$skip) return(chk$state)

  testers <- AlphaSimR::selectInd(
    pop = source$pop,
    nInd = as.integer(cfg$n_testers),
    trait = AlphaSimR::selIndex,
    use = cfg$initial_tester_selection_use,
    b = cfg$initial_tester_selection_weights,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = testers,
    stage = "Female_Testers",
    source = source,
    ready_in_years = 0,
    inherit_genotypes = TRUE,
    selection_strategy = "Fixed female testers for all training and AYT testcrosses",
    stream = "tester"
  )
}

make_training_R_crosses <- function(state, cfg, year) {
  founders <- select_latest_available(state, "R_Founders", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, founders, cfg, event_name = "make_training_R_crosses")
  if (chk$skip) return(chk$state)

  f1 <- make_crosses_from_parents(
    parents = founders$pop,
    n_crosses = cfg$n_training_crosses,
    n_progeny_per_cross = 1L,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = f1,
    stage = "Training_F1",
    source = founders,
    ready_in_years = cfg$training_cross_duration_years,
    cross_strategy = "Elite x elite crosses among the 5 related R founders",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$training_cross_cost_per_cross,
    cost_units = cfg$n_training_crosses,
    cost_event = "crossing",
    cost_unit = "cross"
  )
}

advance_training_RLines_by_SSD <- function(state, cfg, year) {
  f1 <- select_latest_available(state, "Training_F1", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, f1, cfg, event_name = "advance_training_RLines_by_SSD")
  if (chk$skip) return(chk$state)

  rlines <- advance_by_ssd(
    f1_pop = f1$pop,
    n_lines = cfg$n_training_rlines,
    n_selfing_generations = cfg$training_ssd_selfing_generations,
    simParam = state$sim$SP
  )

  state <- put_stage_pop(
    state = state,
    pop = rlines,
    stage = "Training_RLines",
    source = f1,
    ready_in_years = cfg$training_ssd_duration_years,
    selection_strategy = "No selection during SSD/RGA",
    cross_strategy = "SSD/RGA from F1 to inbred R lines",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$training_line_development_cost_per_line,
    cost_event = "line_development",
    cost_unit = "line"
  )
  add_stage_cost(
    state = state,
    event = "seed_increase",
    n_units = nInd(rlines),
    unit_cost = cfg$training_seed_increase_cost_per_line,
    stage = "Training_RLines",
    unit = "line"
  )
}

genotype_training_RLines <- function(state, cfg, year) {
  run_genotyping(
    state,
    list(
      input_stage = "Training_RLines",
      stream = "main",
      input_policy = "latest_one",
      chip = cfg$snpChip,
      duration_years = cfg$genotyping_duration_years,
      cost_per_sample = cfg$genotype_cost_per_sample
    )
  )
}

make_sparse_training_testcrosses <- function(state, cfg, year) {
  rlines <- select_latest_available(state, "Training_RLines", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, rlines, cfg, event_name = "make_sparse_training_testcrosses")
  if (chk$skip) return(chk$state)

  testers <- select_latest_available(state, "Female_Testers", stream = "tester", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, testers, cfg, event_name = "make_sparse_training_testcrosses")
  if (chk$skip) return(chk$state)

  plan <- make_sparse_testcross_plan(
    n_rlines = nInd(rlines$pop),
    n_testers = nInd(testers$pop),
    anchor_fraction = cfg$training_testcross_anchor_fraction
  )
  hybrids <- make_testcross_pop(
    rlines = rlines$pop,
    testers = testers$pop,
    testcross_plan = plan,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = hybrids,
    stage = "Training_Testcrosses",
    source_ids = c(rlines$source_ids, testers$source_ids),
    ready_in_years = cfg$training_testcross_duration_years,
    cross_strategy = "Sparse R-line x tester design with high connectivity across two testers",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$training_testcross_cost_per_cross,
    cost_units = nrow(plan),
    cost_event = "crossing",
    cost_unit = "testcross"
  )
}

run_partner_MET_and_estimate_GCA <- function(state, cfg, year) {
  testcrosses <- select_latest_available(state, "Training_Testcrosses", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, testcrosses, cfg, event_name = "run_partner_MET_and_estimate_GCA")
  if (chk$skip) return(chk$state)

  state <- run_sparse_partner_trial(
    state = state,
    pop = testcrosses$pop,
    output_stage = "Training_Testcross_MET",
    input_cohorts = testcrosses$source_ids,
    cfg = cfg,
    trial_cfg_prefix = "training_met",
    selection_strategy = "Sparse partner MET for GCA, SCA, and GxE estimation",
    stream = "main"
  )

  met <- select_latest_available(
    state = state,
    stage = "Training_Testcross_MET",
    stream = "main",
    n = 1L,
    combine = TRUE,
    include_not_ready = TRUE,
    silent = TRUE
  )
  rlines <- select_latest_available(state, "Training_RLines", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  summaries <- attach_gca_sca_summaries(rlines$pop, met$pop, cfg)
  state$pops[[met$source_ids]] <- summaries$hybrids

  put_stage_pop(
    state = state,
    pop = summaries$rlines,
    stage = "Training_RLines_GCA",
    source = rlines,
    ready_in_years = cfg$training_met_duration_years,
    inherit_genotypes = TRUE,
    selection_strategy = paste(
      "Attach GCA, SCA, and desired-gain index summaries from sparse partner MET",
      paste(met$source_ids, collapse = ";")
    ),
    stream = "main"
  )
}

select_training_R_parents_by_GCA_index <- function(state, cfg, year) {
  rlines <- select_latest_available(state, "Training_RLines_GCA", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, rlines, cfg, event_name = "select_training_R_parents_by_GCA_index")
  if (chk$skip) return(chk$state)

  selected <- bp_select_synthetic(
    pop = rlines$pop,
    n_select = cfg$n_recurrent_parents,
    synthetic_trait = cfg$synthetic_trait,
    use = "pheno",
    state = state,
    varE = cfg$varE,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = selected,
    stage = "R_Parents",
    source = rlines,
    ready_in_years = 0,
    cycle_id = "cycle_0",
    inherit_genotypes = TRUE,
    selection_strategy = "Best 5 R lines by GCA desired-gain index",
    stream = "main"
  )
}

train_fixed_GS_model_from_GCA_index <- function(state, cfg, year) {
  run_train_gp_model(
    state,
    list(
      from_stage = "Training_RLines_GCA",
      stream = "main",
      chip = cfg$snpChip,
      response = "synthetic_pheno",
      synthetic_trait = cfg$synthetic_trait_name,
      lookback_years = cfg$training_model_lookback_years,
      model_id = cfg$gs_model_id
    )
  )
}

make_recurrent_R_crosses <- function(state, cfg, year, cycle) {
  parents <- select_latest_available(state, "R_Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, parents, cfg, event_name = "make_recurrent_R_crosses")
  if (chk$skip) return(chk$state)

  f1 <- make_crosses_from_parents(
    parents = parents$pop,
    n_crosses = cfg$n_recurrent_crosses,
    n_progeny_per_cross = 1L,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = f1,
    stage = "R_F1",
    source = parents,
    ready_in_years = cfg$recurrent_cross_duration_years,
    cycle_id = paste0("cycle_", as.integer(cycle)),
    cross_strategy = "Elite R x elite R crosses among current 5 selected parents",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$recurrent_cross_cost_per_cross,
    cost_units = cfg$n_recurrent_crosses,
    cost_event = "crossing",
    cost_unit = "cross"
  )
}

advance_recurrent_to_F5 <- function(state, cfg, year, cycle) {
  f1 <- select_latest_available(state, "R_F1", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, f1, cfg, event_name = "advance_recurrent_to_F5")
  if (chk$skip) return(chk$state)

  candidates <- advance_by_ssd(
    f1_pop = f1$pop,
    n_lines = cfg$n_recurrent_f5_candidates,
    n_selfing_generations = cfg$recurrent_ssd_selfing_generations,
    simParam = state$sim$SP
  )

  state <- put_stage_pop(
    state = state,
    pop = candidates,
    stage = "R_Candidates_F5",
    source = f1,
    ready_in_years = cfg$recurrent_ssd_duration_years,
    cycle_id = paste0("cycle_", as.integer(cycle)),
    selection_strategy = "No phenotypic selection during recurrent SSD/RGA",
    cross_strategy = "SSD/RGA from recurrent F1 to F5",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$recurrent_line_development_cost_per_line,
    cost_event = "line_development",
    cost_unit = "line"
  )
  add_stage_cost(
    state = state,
    event = "seed_increase",
    n_units = nInd(candidates),
    unit_cost = cfg$recurrent_seed_increase_cost_per_line,
    stage = "R_Candidates_F5",
    unit = "line"
  )
}

genotype_recurrent_F5_candidates <- function(state, cfg, year, cycle) {
  run_genotyping(
    state,
    list(
      input_stage = "R_Candidates_F5",
      stream = "main",
      input_policy = "latest_one",
      include_not_ready = FALSE,
      chip = cfg$snpChip,
      duration_years = cfg$genotyping_duration_years,
      cost_per_sample = cfg$genotype_cost_per_sample
    )
  )
}

predict_and_select_recurrent_R_parents <- function(state, cfg, year, cycle) {
  candidates <- select_latest_available(state, "R_Candidates_F5", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, candidates, cfg, event_name = "predict_and_select_recurrent_R_parents")
  if (chk$skip) return(chk$state)

  scored <- score_latest_candidates_with_fixed_model(state, candidates, cfg, "R_Candidates_F5")
  selected <- bp_select_synthetic(
    pop = scored$scored$pop,
    n_select = cfg$n_recurrent_parents,
    synthetic_trait = cfg$synthetic_trait,
    use = "ebv",
    state = scored$state,
    prediction = "direct",
    simParam = scored$state$sim$SP
  )

  put_stage_pop(
    state = scored$state,
    pop = selected,
    stage = "R_Parents",
    source = candidates,
    ready_in_years = 0,
    cycle_id = paste0("cycle_", as.integer(cycle)),
    inherit_genotypes = TRUE,
    selection_strategy = "Best 5 F5 candidates by fixed-model index GEBV",
    stream = "main"
  )
}

predict_and_select_final30_RLines <- function(state, cfg, year, cycle = 5L) {
  candidates <- select_latest_available(state, "R_Candidates_F5", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, candidates, cfg, event_name = "predict_and_select_final30_RLines")
  if (chk$skip) return(chk$state)

  scored <- score_latest_candidates_with_fixed_model(state, candidates, cfg, "R_Candidates_F5")
  selected <- bp_select_synthetic(
    pop = scored$scored$pop,
    n_select = cfg$n_final_r_lines,
    synthetic_trait = cfg$synthetic_trait,
    use = "ebv",
    state = scored$state,
    prediction = "direct",
    simParam = scored$state$sim$SP
  )

  put_stage_pop(
    state = scored$state,
    pop = selected,
    stage = "R_Final30_F5",
    source = candidates,
    ready_in_years = 0,
    cycle_id = paste0("cycle_", as.integer(cycle)),
    inherit_genotypes = TRUE,
    selection_strategy = sprintf(
      "Best %d F5 candidates by fixed-model index GEBV",
      as.integer(cfg$n_final_r_lines)
    ),
    stream = "main"
  )
}

advance_final30_F5_to_F7 <- function(state, cfg, year) {
  final_f5 <- select_latest_available(state, "R_Final30_F5", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, final_f5, cfg, event_name = "advance_final30_F5_to_F7")
  if (chk$skip) return(chk$state)

  final_f7 <- final_f5$pop
  for (gen in seq_len(as.integer(cfg$final_f5_to_f7_selfing_generations))) {
    final_f7 <- AlphaSimR::self(final_f7, nProgeny = 1L, simParam = state$sim$SP)
  }

  put_stage_pop(
    state = state,
    pop = final_f7,
    stage = "R_Final30_F7",
    source = final_f5,
    ready_in_years = cfg$final_f5_to_f7_duration_years,
    cycle_id = final_f5$cycle_id,
    inherit_genotypes = TRUE,
    selection_strategy = "Advance selected final R lines from F5 to F7",
    cross_strategy = "SSD/RGA F5 to F7",
    stream = "main",
    cost_per_unit = cfg$final_line_development_cost_per_line,
    cost_event = "line_development",
    cost_unit = "line"
  )
}

make_final_testcross_hybrids <- function(state, cfg, year) {
  rlines <- select_latest_available(state, "R_Final30_F7", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, rlines, cfg, event_name = "make_final_testcross_hybrids")
  if (chk$skip) return(chk$state)

  testers <- select_latest_available(state, "Female_Testers", stream = "tester", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, testers, cfg, event_name = "make_final_testcross_hybrids")
  if (chk$skip) return(chk$state)

  plan <- make_full_testcross_plan(nInd(rlines$pop), nInd(testers$pop))
  hybrids <- make_testcross_pop(
    rlines = rlines$pop,
    testers = testers$pop,
    testcross_plan = plan,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = hybrids,
    stage = "AYT_Hybrids",
    source_ids = c(rlines$source_ids, testers$source_ids),
    ready_in_years = cfg$final_testcross_duration_years,
    cycle_id = "cycle_5",
    cross_strategy = sprintf(
      "All %d final R lines crossed to both fixed female testers",
      nInd(rlines$pop)
    ),
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$final_testcross_cost_per_cross,
    cost_units = nrow(plan),
    cost_event = "crossing",
    cost_unit = "testcross"
  )
}

run_partner_AYT <- function(state, cfg, year) {
  hybrids <- select_latest_available(state, "AYT_Hybrids", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, hybrids, cfg, event_name = "run_partner_AYT")
  if (chk$skip) return(chk$state)

  run_phenotype_trial(
    state = state,
    pop = hybrids$pop,
    output_stage = "AYT_Results",
    input_cohorts = hybrids$source_ids,
    selection_strategy = "Partner AYT on all final R x tester hybrid combinations",
    traits = trait_indices(cfg),
    synthetic_traits = cfg$synthetic_trait_name,
    n_loc = cfg$ayt_n_locs,
    reps = cfg$ayt_reps,
    varE = cfg$varE,
    duration_years = cfg$ayt_duration_years,
    stream = "main",
    cycle_id = hybrids$cycle_id,
    cost_per_plot = cfg$ayt_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$ayt_env_means,
    env_mean_sd = cfg$ayt_env_mean_sd,
    env_year_sd = cfg$ayt_env_year_sd,
    log_per_environment = TRUE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

release_hybrids_from_AYT <- function(state, cfg, year) {
  ayt <- select_latest_available(state, "AYT_Results", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, ayt, cfg, event_name = "release_hybrids_from_AYT")
  if (chk$skip) return(chk$state)

  yield_idx <- match("yield", cfg$trait_names)
  if (is.na(yield_idx)) stop("cfg$trait_names must include 'yield'.", call. = FALSE)
  keep <- which(ayt$pop@pheno[, yield_idx] > cfg$ayt_local_check_yield)

  if (!length(keep)) {
    return(bp_log_event(
      state = state,
      fn = "release_hybrids_from_AYT",
      event_type = "release",
      stage = "Released_Hybrids",
      source_ids = ayt$source_ids,
      output_id = NA_character_,
      event_string = "No AYT hybrids outperformed the local commercial yield check.",
      template_string = "Release hybrids above local check",
      details = list(local_check_yield = cfg$ayt_local_check_yield, n_released = 0L)
    ))
  }

  put_stage_pop(
    state = state,
    pop = ayt$pop[keep],
    stage = "Released_Hybrids",
    source = ayt,
    ready_in_years = 0,
    cycle_id = ayt$cycle_id,
    inherit_genotypes = TRUE,
    selection_strategy = "Observed AYT yield strictly outperforms local commercial check",
    stream = "main"
  )
}

# Runners ------------------------------------------------------------------

run_simulation <- function(state, results_base = data.frame(), cfg, results_file = NULL, seed = NULL, start_year = NULL) {
  if (!is.null(seed)) set.seed(seed)

  report_start_year <- if (is.null(start_year)) state$time$t else start_year
  if (is.null(state$sim$synthetic_traits[[cfg$synthetic_trait_name]])) {
    state <- bp_register_synthetic_traits(state, cfg$synthetic_trait)
  }

  results <- record_yearly_outputs(
    results = data.frame(),
    results_base = results_base,
    state = state,
    year = state$time$t - report_start_year,
    cfg = cfg
  )

  state <- initialize_R_founders(state, cfg, state$time$t)
  state <- initialize_fixed_testers(state, cfg, state$time$t)

  state <- make_training_R_crosses(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$training_cross_duration_years)

  state <- advance_training_RLines_by_SSD(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$training_ssd_duration_years)

  state <- genotype_training_RLines(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$genotyping_duration_years)

  state <- make_sparse_training_testcrosses(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$training_testcross_duration_years)

  state <- run_partner_MET_and_estimate_GCA(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$training_met_duration_years)

  state <- select_training_R_parents_by_GCA_index(state, cfg, state$time$t)
  state <- train_fixed_GS_model_from_GCA_index(state, cfg, state$time$t)

  results <- record_yearly_outputs(
    results = results,
    results_base = results_base,
    state = state,
    year = state$time$t - report_start_year,
    cfg = cfg
  )

  for (cycle in seq_len(as.integer(cfg$n_recurrent_cycles - 1L))) {
    state <- make_recurrent_R_crosses(state, cfg, state$time$t, cycle)
    state <- bp_advance_time_years(state, cfg$recurrent_cross_duration_years)

    state <- advance_recurrent_to_F5(state, cfg, state$time$t, cycle)
    state <- bp_advance_time_years(state, cfg$recurrent_ssd_duration_years)

    state <- genotype_recurrent_F5_candidates(state, cfg, state$time$t, cycle)
    state <- bp_advance_time_years(state, cfg$genotyping_duration_years)

    state <- predict_and_select_recurrent_R_parents(state, cfg, state$time$t, cycle)
    results <- record_yearly_outputs(
      results = results,
      results_base = results_base,
      state = state,
      year = state$time$t - report_start_year,
      cfg = cfg
    )
  }

  final_cycle <- as.integer(cfg$n_recurrent_cycles)
  state <- make_recurrent_R_crosses(state, cfg, state$time$t, final_cycle)
  state <- bp_advance_time_years(state, cfg$recurrent_cross_duration_years)

  state <- advance_recurrent_to_F5(state, cfg, state$time$t, final_cycle)
  state <- bp_advance_time_years(state, cfg$recurrent_ssd_duration_years)

  state <- genotype_recurrent_F5_candidates(state, cfg, state$time$t, final_cycle)
  state <- bp_advance_time_years(state, cfg$genotyping_duration_years)

  state <- predict_and_select_final30_RLines(state, cfg, state$time$t, final_cycle)
  state <- advance_final30_F5_to_F7(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$final_f5_to_f7_duration_years)

  state <- make_final_testcross_hybrids(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$final_testcross_duration_years)

  state <- run_partner_AYT(state, cfg, state$time$t)
  state <- bp_advance_time_years(state, cfg$ayt_duration_years)

  state <- release_hybrids_from_AYT(state, cfg, state$time$t)
  results <- record_yearly_outputs(
    results = results,
    results_base = results_base,
    state = state,
    year = state$time$t - report_start_year,
    cfg = cfg
  )

  if (!is.null(results_file)) write.csv(results, file = results_file, row.names = FALSE)
  invisible(list(state = state, results = results))
}
