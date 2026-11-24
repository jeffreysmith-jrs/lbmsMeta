#################################################
## 05_summarise_models.R
##
## Summarises the fitted models (steps 03-04) into
## tables used by step 06 and the manuscript text.
## Every summary is on the log response ratio scale;
## fold_change = exp(median), percent_change =
## 100 * (fold_change - 1). Moderate support: the 80%
## credible interval excludes zero; strong support:
## the 95% interval excludes zero.
##
## Output: output/summaries/
##   all_model_coefficients.csv  every fixed effect of every model
##   model_diagnostics.csv       R-hat and ESS of every model
##   main_effects.csv            main-model intercepts
##   all_<moderator>_effects.csv estimates per moderator level
##   all_<analysis>_effects.csv  sensitivity-analysis estimates
##   all_area_curves.csv, all_time_curves.csv
##   all_study_heterogeneity.csv, all_tau.csv,
##   all_leave_one_out.csv, all_model_counts.csv
#################################################

library(brms)
library(posterior)
library(dplyr)

paths_file <- Find(file.exists, c("R/paths.R", file.path(Sys.getenv("LBMS_ANALYSIS_DIR"), "R", "paths.R")))
if(is.null(paths_file)) stop("Run from the repository root (the folder that contains run_analysis.R).")
source(paths_file)
source(file.path(code_dir, "R", "helpers.R"))
source(file.path(code_dir, "R", "models.R"))

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

## Models summarised, by folder under
## output/models/<intervention>/
moderator_models <- c(
  "main_model", "taxa_model", "zone_model", "biome_model", "structure_model",
  "human_use_model", "experiment_type_model", "area_model", "time_model",
  "forestation_type_model"
)

sensitivity_models <- c(
  "mean_aggregation_model", "aggregated_aggregation_model",
  "weighted_model", "unweighted_model",
  "pet_model", "peese_model",
  "excl_high_n_model", "one_per_study_model",
  "prior_tight_model", "prior_wide_model"
)

## Saved model files for one intervention, in
## alphabetical order within each folder
model_files_for <- function(intervention){
  dir <- file.path(results_dir, intervention)
  top <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  sens <- list.files(file.path(dir, "sensitivity"), pattern = "\\.rds$", full.names = TRUE)
  c(top[model_name_of(top) %in% moderator_models],
    sens[model_name_of(sens) %in% sensitivity_models])
}

model_name_of <- function(f) tools::file_path_sans_ext(basename(f))

#################################################
## Posterior summaries
#################################################

extract_parameter_summary <- function(draws, parameter){

  vals <- draws[[parameter]]

  data.frame(
    parameter = parameter,
    mean = mean(vals, na.rm = TRUE),
    median = median(vals, na.rm = TRUE),
    sd = sd(vals, na.rm = TRUE),
    q025 = quantile(vals, .025),
    q05 = quantile(vals, .05),
    q10 = quantile(vals, .10),
    q90 = quantile(vals, .90),
    q95 = quantile(vals, .95),
    q975 = quantile(vals, .975)
  )

}

## Every fixed effect (b_*) of one model
process_model <- function(model, intervention, model_name){

  draws <- posterior::as_draws_df(model)
  parameters <- grep("^b_", names(draws), value = TRUE)

  if(length(parameters) == 0) return(NULL)

  out <- lapply(parameters, function(p){
    tmp <- extract_parameter_summary(draws, p)
    tmp$intervention <- intervention
    tmp$model <- model_name
    tmp
  })

  do.call(rbind, out)

}

summarise_effect_draws <- function(draws, intervention, level){

  data.frame(
    intervention = intervention,
    level = level,
    median = median(draws),
    q025 = quantile(draws, .025),
    q05 = quantile(draws, .05),
    q10 = quantile(draws, .10),
    q90 = quantile(draws, .90),
    q95 = quantile(draws, .95),
    q975 = quantile(draws, .975)
  )

}

add_effect_metrics <- function(df){

  df$fold_change <- exp(df$median)
  df$percent_change <- (exp(df$median) - 1) * 100
  df$moderate_support <- sign(df$q10) == sign(df$q90)
  df$strong_support <- sign(df$q025) == sign(df$q975)
  df

}

## TRUE if the 95% credible interval excludes zero
strong_support <- function(draws){
  sign(quantile(draws, 0.025, na.rm = TRUE)) == sign(quantile(draws, 0.975, na.rm = TRUE))
}

#################################################
## Convergence diagnostics
##
## The rule used while fitting (R/models.R): max
## R-hat <= 1.01, min bulk and tail ESS >= 400, no
## divergent transitions.
#################################################

extract_diagnostics <- function(model){
  d <- model_convergence_diagnostics(model)
  data.frame(max_rhat = d$max_rhat, min_bulk_ess = d$min_ess_bulk, min_tail_ess = d$min_ess_tail,
             n_divergent = d$n_divergent, converged = d$converged)
}

#################################################
## Study-level heterogeneity (main model)
##
## effects: the pooled intercept ("Overall") and each
## study's partially pooled effect (intercept + that
## study's random intercept), with the ratio of each
## study's effect to the pooled effect,
## exp(r_study_id).
##
## tau: the between-study SD (tau); I2 = tau^2 /
## (tau^2 + sigma^2), the share of residual variance
## between studies (Nakagawa & Santos 2012; the main
## model is unweighted, so sigma^2 is the residual
## variance, not the sampling variance).
#################################################

extract_study_heterogeneity <- function(model, intervention){

  draws <- posterior::as_draws_df(model)

  re_cols <- grep("^r_study_id\\[.*,Intercept\\]$", names(draws), value = TRUE)

  if(!"b_Intercept" %in% names(draws) || length(re_cols) == 0) return(NULL)

  overall_row <- extract_parameter_summary(draws, "b_Intercept")
  overall_row$study_id <- "Overall"
  overall_row$row_type <- "Overall"

  study_df <- do.call(rbind, lapply(re_cols, function(col){

    study_draws <- draws$b_Intercept + draws[[col]]
    ratio_draws <- exp(draws[[col]])

    data.frame(
      parameter = col,
      mean = mean(study_draws, na.rm = TRUE),
      median = median(study_draws, na.rm = TRUE),
      sd = sd(study_draws, na.rm = TRUE),
      q025 = quantile(study_draws, .025),
      q05 = quantile(study_draws, .05),
      q10 = quantile(study_draws, .10),
      q90 = quantile(study_draws, .90),
      q95 = quantile(study_draws, .95),
      q975 = quantile(study_draws, .975),
      study_id = sub("^r_study_id\\[(.*),Intercept\\]$", "\\1", col),
      row_type = "Study",
      ratio_to_fixed_effect_median = median(ratio_draws, na.rm = TRUE),
      ratio_to_fixed_effect_q025 = quantile(ratio_draws, .025),
      ratio_to_fixed_effect_q975 = quantile(ratio_draws, .975)
    )

  }))

  overall_row$ratio_to_fixed_effect_median <- NA_real_
  overall_row$ratio_to_fixed_effect_q025 <- NA_real_
  overall_row$ratio_to_fixed_effect_q975 <- NA_real_

  effects_df <- rbind(overall_row[, names(study_df)], study_df)
  effects_df$intervention <- intervention

  tau_draws <- draws$sd_study_id__Intercept

  tau_df <- data.frame(
    intervention = intervention,
    tau_median = median(tau_draws, na.rm = TRUE),
    tau_mean = mean(tau_draws, na.rm = TRUE),
    tau_q025 = quantile(tau_draws, .025),
    tau_q975 = quantile(tau_draws, .975),
    n_studies = length(re_cols),
    n_observations = nrow(model$data)
  )

  i2_draws <- tau_draws^2 / (tau_draws^2 + draws$sigma^2)

  tau_df$i2_median <- median(i2_draws, na.rm = TRUE) * 100
  tau_df$i2_q025 <- quantile(i2_draws, .025) * 100
  tau_df$i2_q975 <- quantile(i2_draws, .975) * 100

  list(effects = effects_df, tau = tau_df)

}

#################################################
## Area and age response curves
##
## Intercept + slope * x across the intervention's
## own observed range of log(area) or age (the
## S_aggregation term is at its average, zero).
#################################################

build_prediction_curve <- function(model, predictor_grid, predictor_prefix){

  draws <- posterior::as_draws_df(model)

  slope_name <- grep(paste0("^b_", predictor_prefix), names(draws), value = TRUE)

  if(!"b_Intercept" %in% names(draws) || length(slope_name) == 0) return(NULL)

  intercept <- draws$b_Intercept
  slope <- draws[[slope_name[1]]]

  do.call(rbind, lapply(predictor_grid, function(x){
    pred <- intercept + slope * x
    data.frame(
      predictor = x,
      median = median(pred),
      lower80 = quantile(pred, .10),
      upper80 = quantile(pred, .90),
      lower95 = quantile(pred, .025),
      upper95 = quantile(pred, .975)
    )
  }))

}

observed_grid <- function(x){
  seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
}

#################################################
## Coefficients, diagnostics, heterogeneity, curves
#################################################

all_coefficients <- data.frame()
all_diagnostics <- data.frame()
all_area_curves <- data.frame()
all_time_curves <- data.frame()
all_study_heterogeneity <- data.frame()
all_tau <- data.frame()

for(intervention in interventions){

  for(f in model_files_for(intervention)){

    model_name <- model_name_of(f)
    message(intervention, " : ", model_name)

    model <- tryCatch(readRDS(f), error = function(e) NULL)
    if(is.null(model)) next

    tmp <- tryCatch(process_model(model, intervention, model_name), error = function(e) NULL)
    if(!is.null(tmp)) all_coefficients <- rbind(all_coefficients, tmp)

    diag_tmp <- tryCatch(extract_diagnostics(model), error = function(e) NULL)

    if(!is.null(diag_tmp)){

      diag_tmp$intervention <- intervention
      diag_tmp$model <- model_name

      all_diagnostics <- rbind(all_diagnostics, diag_tmp)

    }

    if(model_name == "main_model"){

      het <- tryCatch(extract_study_heterogeneity(model, intervention), error = function(e) NULL)

      if(!is.null(het)){
        all_study_heterogeneity <- rbind(all_study_heterogeneity, het$effects)
        all_tau <- rbind(all_tau, het$tau)
      }

    }

    if(model_name == "area_model"){
      curve <- build_prediction_curve(model, observed_grid(log(model$data$mean_size_ha)), "logmean_size_ha")
      if(!is.null(curve)){
        curve$intervention <- intervention
        all_area_curves <- rbind(all_area_curves, curve)
      }
    }

    if(model_name == "time_model"){
      curve <- build_prediction_curve(model, observed_grid(model$data$intervention_length), "intervention_length")
      if(!is.null(curve)){
        curve$intervention <- intervention
        all_time_curves <- rbind(all_time_curves, curve)
      }
    }

  }

}

## Pooled forestation models: diagnostics only (their
## estimates are summarised per moderator below)
for(f in model_files_for("forestation")){
  diag_tmp <- tryCatch(extract_diagnostics(readRDS(f)), error = function(e) NULL)
  if(!is.null(diag_tmp)){
    diag_tmp$intervention <- "forestation"
    diag_tmp$model <- model_name_of(f)
    all_diagnostics <- rbind(all_diagnostics, diag_tmp)
  }
}

for(i in which(!all_diagnostics$converged)){
  log_message("NOT CONVERGED: ", all_diagnostics$intervention[i], " / ", all_diagnostics$model[i],
              " (max R-hat ", round(all_diagnostics$max_rhat[i], 4), ", min ESS bulk/tail ",
              round(all_diagnostics$min_bulk_ess[i]), "/", round(all_diagnostics$min_tail_ess[i]),
              ", divergent ", all_diagnostics$n_divergent[i], ")")
}

all_coefficients <- add_effect_metrics(all_coefficients)

#################################################
## Leave-one-study-out (written by step 04)
#################################################

all_leave_one_out <- bind_rows(lapply(interventions, function(intervention){
  f <- file.path(results_dir, intervention, "sensitivity", "leave_one_out_summary.csv")
  if(file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}))

if(nrow(all_leave_one_out) > 0){

  all_leave_one_out$fold_change <- exp(all_leave_one_out$median)
  all_leave_one_out$percent_change <- (exp(all_leave_one_out$median) - 1) * 100

  write.csv(all_leave_one_out, file.path(summary_dir, "all_leave_one_out.csv"), row.names = FALSE)

}

#################################################
## Estimates per moderator level
##
## Each level's estimate is the intercept (reference
## level) or intercept + its coefficient; the
## reference level is the first level of the factor
## in the model's data. Two-level moderators also get
## a "Difference" row (the coefficient itself).
## Multi-level moderators get compact letters
## instead: two levels share a letter unless the
## 95% CI of their difference excludes zero (strong
## support).
#################################################

## Compact-letter display for a named list of draws.
## Groups are taken in order of decreasing median, so
## the highest gets "A".
assign_cld_letters <- function(draws_list){

  groups <- names(draws_list)
  n <- length(groups)

  diff_mat <- matrix(FALSE, n, n, dimnames = list(groups, groups))

  for(i in seq_len(n)){
    for(j in seq_len(n)){
      if(i != j) diff_mat[i, j] <- strong_support(draws_list[[i]] - draws_list[[j]])
    }
  }

  groups_ord <- groups[order(sapply(draws_list, median), decreasing = TRUE)]

  letter_groups <- list()

  for(g in groups_ord){

    placed <- FALSE

    ## g joins every existing group whose members it is
    ## not strongly different from
    for(k in seq_along(letter_groups)){
      members <- letter_groups[[k]]
      if(all(!diff_mat[g, members])){
        letter_groups[[k]] <- c(members, g)
        placed <- TRUE
      }
    }

    if(!placed) letter_groups[[length(letter_groups) + 1]] <- g

  }

  labels <- setNames(rep("", n), groups)

  for(k in seq_along(letter_groups)){
    for(g in letter_groups[[k]]){
      labels[g] <- paste0(labels[g], LETTERS[k])
    }
  }

  labels

}

## Posterior draws of each level, reference first
level_draws_for <- function(model, variable){

  draws <- posterior::as_draws_df(model)
  if(!"b_Intercept" %in% names(draws)) return(NULL)

  ## Levels are matched to coefficients by position,
  ## because brms drops spaces from coefficient names
  levels_present <- levels(factor(model$data[[variable]]))
  pars <- grep(paste0("^b_", variable), names(draws), value = TRUE)

  out <- list()
  out[[levels_present[1]]] <- draws$b_Intercept
  for(i in seq_along(pars)){
    out[[levels_present[i + 1]]] <- draws$b_Intercept + draws[[pars[i]]]
  }

  attr(out, "coefficients") <- lapply(pars, function(p) draws[[p]])
  out

}

two_level_effects <- function(model_name, variable, intervention_set){

  out <- data.frame()

  for(intervention in intervention_set){

    f <- file.path(results_dir, intervention, paste0(model_name, ".rds"))
    if(!file.exists(f)) next

    lv <- level_draws_for(readRDS(f), variable)
    if(is.null(lv)) next

    for(level in names(lv)){
      out <- bind_rows(out, summarise_effect_draws(lv[[level]], intervention, level))
    }

    coefs <- attr(lv, "coefficients")
    if(length(coefs) > 0){
      out <- bind_rows(out, summarise_effect_draws(coefs[[1]], intervention, "Difference"))
    }

  }

  add_effect_metrics(out)

}

multi_level_effects <- function(model_name, variable, intervention_set){

  out <- data.frame()

  for(intervention in intervention_set){

    f <- file.path(results_dir, intervention, paste0(model_name, ".rds"))
    if(!file.exists(f)) next

    lv <- level_draws_for(readRDS(f), variable)
    if(is.null(lv)) next

    letters_by_level <-
      if(length(lv) >= 2) assign_cld_letters(lv) else setNames(rep("", length(lv)), names(lv))

    for(level in names(lv)){
      out <- bind_rows(
        out,
        summarise_effect_draws(lv[[level]], intervention, level) %>%
          mutate(is_reference = (level == names(lv)[1]), cld_letter = letters_by_level[[level]])
      )
    }

  }

  add_effect_metrics(out)

}

with_forestation <- c(interventions, "forestation")

taxa_effects <- multi_level_effects("taxa_model", "taxa_grouped", interventions)
taxa_effects$level <- paste0(toupper(substring(taxa_effects$level, 1, 1)), substring(taxa_effects$level, 2))

biome_effects <- two_level_effects("biome_model", "forested", with_forestation)
zone_effects <- two_level_effects("zone_model", "zone", with_forestation)
structure_effects <- two_level_effects("structure_model", "reference_habitat_structure_axis", with_forestation)
human_use_effects <- multi_level_effects("human_use_model", "reference_habitat_human_use_axis", with_forestation)
experiment_type_effects <- two_level_effects("experiment_type_model", "experiment_type", with_forestation)

#################################################
## Forestation method x biome
##
## Cell means of the forestation_type_model
## interaction, with compact letters across the four
## cells. "natural forest restoration" is reported
## as Natural forest regrowth and "timber
## plantation" as Plantation.
#################################################

forestation_method_effects <- data.frame()

forestation_type_model_file <- file.path(results_dir, "forestation", "forestation_type_model.rds")

if(file.exists(forestation_type_model_file)){

  draws <- posterior::as_draws_df(readRDS(forestation_type_model_file))
  b_pars <- grep("^b_", names(draws), value = TRUE)

  forest_type_par <- b_pars[grepl("^b_forest_type", b_pars) & !grepl(":", b_pars)]
  forested_par <- b_pars[grepl("^b_forested", b_pars) & !grepl(":", b_pars)]
  interaction_par <- b_pars[grepl(":", b_pars)]

  cell_draws <- list(
    "Natural forest regrowth.Forested" = draws$b_Intercept,
    "Natural forest regrowth.Open" = draws$b_Intercept + draws[[forested_par]],
    "Plantation.Forested" = draws$b_Intercept + draws[[forest_type_par]],
    "Plantation.Open" =
      draws$b_Intercept + draws[[forested_par]] + draws[[forest_type_par]] + draws[[interaction_par]]
  )

  cld_letters <- assign_cld_letters(cell_draws)

  for(cell in names(cell_draws)){
    parts <- strsplit(cell, "\\.")[[1]]
    row <- summarise_effect_draws(cell_draws[[cell]], "forestation", parts[2])
    row$method <- parts[1]
    row$cld_letter <- cld_letters[cell]
    forestation_method_effects <- bind_rows(forestation_method_effects, row)
  }

  forestation_method_effects <- add_effect_metrics(forestation_method_effects)

}

#################################################
## Sample sizes behind every estimate
##
## Observations and studies each model was fitted on
## (after its own filters), overall, per level of each
## categorical predictor, and per combination when a
## model has two. Used for the "N = observations
## (studies)" labels in the figures.
#################################################

count_model_data <- function(model, intervention, model_name){

  d <- model$data
  if(is.null(d) || !"study_id" %in% names(d)) return(NULL)

  count_rows <- function(df, variable, level){
    data.frame(
      intervention = intervention, model = model_name,
      variable = variable, level = level,
      n_observations = nrow(df),
      n_studies = length(unique(df$study_id)),
      stringsAsFactors = FALSE
    )
  }

  response <- all.vars(stats::formula(model)$formula)[1]
  categorical <- names(d)[vapply(d, function(x) is.factor(x) || is.character(x), logical(1))]
  categorical <- setdiff(categorical, c("study_id", "S_aggregation", response))

  out <- list(count_rows(d, "overall", "overall"))

  for(v in categorical){
    for(lv in sort(unique(as.character(d[[v]])))){
      out[[length(out) + 1]] <- count_rows(d[as.character(d[[v]]) == lv, , drop = FALSE], v, lv)
    }
  }

  if(length(categorical) == 2){
    combo <- paste(as.character(d[[categorical[1]]]), as.character(d[[categorical[2]]]), sep = " | ")
    for(lv in sort(unique(combo))){
      out[[length(out) + 1]] <-
        count_rows(d[combo == lv, , drop = FALSE], paste(categorical, collapse = " | "), lv)
    }
  }

  do.call(rbind, out)

}

model_counts <- list()

for(folder in intersect(list.dirs(results_dir, recursive = FALSE, full.names = FALSE), with_forestation)){
  for(f in model_files_for(folder)){
    m <- tryCatch(readRDS(f), error = function(e) NULL)
    if(is.null(m) || !inherits(m, "brmsfit")) next
    model_counts[[length(model_counts) + 1]] <-
      tryCatch(count_model_data(m, folder, model_name_of(f)), error = function(e) NULL)
  }
}

#################################################
## Write
#################################################

write_summary <- function(df, name){
  write.csv(df, file.path(summary_dir, paste0(name, ".csv")), row.names = FALSE)
}

coefficients_for <- function(models) all_coefficients %>% filter(model %in% models)

write_summary(all_coefficients, "all_model_coefficients")
write_summary(all_diagnostics, "model_diagnostics")
write_summary(coefficients_for("main_model") %>% filter(parameter == "b_Intercept"), "main_effects")

write_summary(taxa_effects, "all_taxa_effects")
write_summary(zone_effects, "all_zone_effects")
write_summary(biome_effects, "all_biome_effects")
write_summary(structure_effects, "all_structure_effects")
write_summary(human_use_effects, "all_human_use_effects")
write_summary(experiment_type_effects, "all_experiment_type_effects")
write_summary(forestation_method_effects, "all_forestation_method_effects")

write_summary(coefficients_for(c("mean_aggregation_model", "aggregated_aggregation_model")), "all_aggregation_effects")
write_summary(coefficients_for("weighted_model"), "all_weighted_effects")
write_summary(coefficients_for("unweighted_model"), "all_unweighted_effects")
write_summary(coefficients_for("pet_model"), "all_pet_effects")
write_summary(coefficients_for("peese_model"), "all_peese_effects")
write_summary(coefficients_for("excl_high_n_model"), "all_high_n_excl_effects")
write_summary(coefficients_for("one_per_study_model"), "all_one_per_study_effects")
write_summary(coefficients_for(c("prior_tight_model", "prior_wide_model")), "all_prior_sensitivity_effects")

if(nrow(all_area_curves) > 0) write_summary(all_area_curves, "all_area_curves")
if(nrow(all_time_curves) > 0) write_summary(all_time_curves, "all_time_curves")

if(nrow(all_study_heterogeneity) > 0){
  write_summary(add_effect_metrics(all_study_heterogeneity), "all_study_heterogeneity")
}
if(nrow(all_tau) > 0) write_summary(all_tau, "all_tau")

write_summary(do.call(rbind, model_counts), "all_model_counts")

log_message(sum(all_diagnostics$converged), " of ", nrow(all_diagnostics),
            " models converged (R-hat <= ", rhat_max, ", ESS >= ", ess_min, ", no divergent transitions).")

log_message("Model summaries complete.")
