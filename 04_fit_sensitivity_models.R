#################################################
## 04_fit_sensitivity_models.R
##
## Sensitivity analyses of the main effect, for each
## intervention:
##   - aggregation method: intercept-only models on
##     the studies reporting mean richness and on
##     those reporting aggregated richness;
##   - variance weighting: on observations with a
##     usable standard error, weighted by it
##     (se(lrr_se_richness)) and unweighted;
##   - publication bias: PET (precision-effect test,
##     slope on the standard error) and PEESE
##     (slope on the variance), both weighted
##     (Stanley & Doucouliagos 2014);
##   - excluding the top 10% of studies by number of
##     observations;
##   - one effect size per study (each study's mean
##     lrr_richness, no random effect);
##   - leave-one-study-out refits of the main model;
##   - prior sensitivity: the main model refitted
##     with tighter and wider priors.
## The leave-one-out and prior refits start from the
## main model fitted in step 03.
##
## Output: output/models/<intervention>/sensitivity/
##   <model>.rds and leave_one_out_summary.csv
#################################################

library(brms)
library(cmdstanr)
library(posterior)

paths_file <- Find(file.exists, c("R/paths.R", file.path(Sys.getenv("LBMS_ANALYSIS_DIR"), "R", "paths.R")))
if(is.null(paths_file)) stop("Run from the repository root (the folder that contains run_analysis.R).")
source(paths_file)
source(file.path(code_dir, "R", "helpers.R"))
source(file.path(code_dir, "R", "models.R"))

## Log row counts and check convergence for every fit
brm <- brm_logged

## The rows of the main model (step 03): observations
## with an effect size in modelled_rows() (R/models.R),
## so Chao, rarefied and other estimates, and summed
## richness at unequal sampling effort, are left out of
## every sensitivity analysis as well. S_aggregation
## keeps its raw spelling here; each model simplifies
## it as needed.
read_effect_data <- function(intervention){
  df <- read.csv(file.path(effects_dir, paste0(intervention, ".csv")), stringsAsFactors = FALSE)
  df <- df[modelled_rows(df), ]
  df[!is.na(df$lrr_richness), ]
}

## Shared brm() arguments. The log rows are labelled
## with the fit_*_model() function that called this.
fit_brm <- function(formula, df, prior){

  brm(
    formula,
    data = df,
    family = gaussian(),
    prior = prior,
    chains = n_chains,
    cores = n_cores,
    iter = n_iter,
    warmup = n_warmup,
    control = brms_control_default,
    backend = "cmdstanr",
    seed = project_seed,
    caller = deparse(sys.call(-1)[[1]])
  )

}

## Observations with a usable standard error (studies
## reporting mean richness with an uncertainty term)
with_se <- function(df){
  df[!is.na(df$lrr_se_richness) & df$lrr_se_richness > 0, ]
}

#################################################
## Aggregation method
#################################################

aggregation_subset <- function(df, method){
  df$S_aggregation <- simplify_aggregation(df$S_aggregation)
  df[!is.na(df$S_aggregation) & df$S_aggregation == method, ]
}

fit_mean_aggregation_model <- function(df){
  df <- aggregation_subset(df, "mean")
  if(nrow(df) < 10) return(NULL)
  fit_brm(lrr_richness ~ 1 + (1 | study_id), df, priors_intercept)
}

fit_aggregated_aggregation_model <- function(df){
  df <- aggregation_subset(df, "aggregated")
  if(nrow(df) < 10) return(NULL)
  fit_brm(lrr_richness ~ 1 + (1 | study_id), df, priors_intercept)
}

#################################################
## Variance weighting
#################################################

fit_weighted_model <- function(df){
  df <- with_se(df)
  if(nrow(df) < 20) return(NULL)
  fit_brm(lrr_richness | se(lrr_se_richness) ~ 1 + (1 | study_id), df, priors_intercept)
}

## Same observations as fit_weighted_model(), not weighted
fit_unweighted_model <- function(df){
  df <- with_se(df)
  if(nrow(df) < 20) return(NULL)
  fit_brm(lrr_richness ~ 1 + (1 | study_id), df, priors_intercept)
}

#################################################
## Publication bias (PET-PEESE)
##
## The slope on the standard error (PET) tests for
## small-study effects (funnel-plot asymmetry). The
## intercepts are the effect extrapolated to a
## perfectly precise study.
#################################################

fit_pet_model <- function(df){
  df <- with_se(df)
  if(nrow(df) < 20) return(NULL)
  fit_brm(lrr_richness | se(lrr_se_richness) ~ 1 + lrr_se_richness + (1 | study_id),
          df, priors_moderator)
}

fit_peese_model <- function(df){
  df <- with_se(df)
  df <- df[!is.na(df$lrr_var_richness), ]
  if(nrow(df) < 20) return(NULL)
  fit_brm(lrr_richness | se(lrr_se_richness) ~ 1 + lrr_var_richness + (1 | study_id),
          df, priors_moderator)
}

#################################################
## Pseudo-replication
#################################################

## Without the top 10% of studies by number of
## observations (high_n_study, step 01)
fit_excl_high_n_model <- function(df){
  df <- df[!is.na(df$high_n_study) & !df$high_n_study, ]
  if(nrow(df) < 10) return(NULL)
  fit_brm(lrr_richness ~ 1 + (1 | study_id), df, priors_intercept)
}

## One row per study: its mean lrr_richness. With one
## row per study there is no study random effect.
fit_one_per_study_model <- function(df){
  study_df <- aggregate(lrr_richness ~ study_id, data = df, FUN = mean)
  if(nrow(study_df) < 10) return(NULL)
  fit_brm(lrr_richness ~ 1, study_df, prior(normal(0, 5), class = "Intercept"))
}

#################################################
## Leave-one-study-out
##
## The main model refitted once per study with that
## study's rows removed. A model read from disk has no
## compiled Stan program, so the first refit compiles
## it and the later refits reuse that program (only
## the data change). Returns one row per excluded
## study with the posterior median and 95% CI of the
## intercept.
#################################################

fit_leave_one_out_models <- function(main_model, intervention){

  base_data <- main_model$data
  study_ids <- unique(base_data$study_id)

  template <- main_model
  compiled <- FALSE

  rows <- lapply(study_ids, function(excluded_study){

    held_out_data <- base_data[base_data$study_id != excluded_study, ]

    if(nrow(held_out_data) < 10) return(NULL)

    refit <- tryCatch(
      update(
        template,
        newdata = held_out_data,
        recompile = !compiled,
        chains = n_chains,
        cores = n_cores,
        iter = n_iter,
        warmup = n_warmup,
        control = brms_control_default,
        seed = project_seed,
        refresh = 0
      ),
      error = function(e){
        log_message("Leave-one-out refit failed for ", intervention, " without ", excluded_study, ": ", conditionMessage(e))
        NULL
      }
    )

    if(is.null(refit)) return(NULL)

    if(!compiled){
      template <<- refit
      compiled <<- TRUE
    }

    draws <- posterior::as_draws_df(refit)

    data.frame(
      intervention = intervention,
      excluded_study_id = excluded_study,
      median = median(draws$b_Intercept, na.rm = TRUE),
      q025 = quantile(draws$b_Intercept, .025),
      q975 = quantile(draws$b_Intercept, .975),
      n_remaining_studies = length(study_ids) - 1,
      n_remaining_observations = nrow(held_out_data)
    )

  })

  out <- do.call(rbind, rows)

  if(is.null(out) || nrow(out) < length(study_ids)){
    log_message("WARNING: leave-one-out for ", intervention, ": ", if(is.null(out)) 0 else nrow(out),
                " of ", length(study_ids), " refits succeeded")
  }

  out

}

#################################################
## Prior sensitivity
##
## The main model (same formula and data) with the
## tight and wide prior schemes (R/models.R). A new
## prior changes the Stan program, so these are fresh
## fits rather than updates.
#################################################

fit_prior_tight_model <- function(main_model){
  fit_brm(formula_main, main_model$data, priors_tight)
}

fit_prior_wide_model <- function(main_model){
  fit_brm(formula_main, main_model$data, priors_wide)
}

#################################################
## Run
#################################################

## A model whose fit failed is NULL: nothing is saved,
## and any file from an earlier run is removed
save_model <- function(model, output_dir, model_name){
  f <- file.path(output_dir, paste0(model_name, ".rds"))
  if(!is.null(model)) saveRDS(model, f) else if(file.exists(f)) file.remove(f)
}

run_sensitivity_models <- function(intervention){

  log_message("Sensitivity: ", intervention)

  df <- read_effect_data(intervention)
  set_model_base(intervention, df)

  output_dir <- file.path(results_dir, intervention, "sensitivity")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  save_model(fit_mean_aggregation_model(df), output_dir, "mean_aggregation_model")
  save_model(fit_aggregated_aggregation_model(df), output_dir, "aggregated_aggregation_model")
  save_model(fit_weighted_model(df), output_dir, "weighted_model")
  save_model(fit_unweighted_model(df), output_dir, "unweighted_model")
  save_model(fit_pet_model(df), output_dir, "pet_model")
  save_model(fit_peese_model(df), output_dir, "peese_model")
  save_model(fit_excl_high_n_model(df), output_dir, "excl_high_n_model")
  save_model(fit_one_per_study_model(df), output_dir, "one_per_study_model")

  ## The refits below need the main model from step 03
  main_model_file <- file.path(results_dir, intervention, "main_model.rds")

  if(!file.exists(main_model_file)){
    log_message(intervention, ": no main_model.rds - leave-one-out and prior refits skipped")
    return(invisible(NULL))
  }

  main_model <- readRDS(main_model_file)

  loo_summary <- tryCatch(
    fit_leave_one_out_models(main_model, intervention),
    error = function(e){
      log_message("Leave-one-out failed for ", intervention, ": ", e$message)
      NULL
    }
  )

  loo_file <- file.path(output_dir, "leave_one_out_summary.csv")
  if(file.exists(loo_file)) file.remove(loo_file)
  if(!is.null(loo_summary) && nrow(loo_summary) > 0){
    write.csv(loo_summary, loo_file, row.names = FALSE)
  }

  save_model(tryCatch(fit_prior_tight_model(main_model), error = function(e) NULL),
             output_dir, "prior_tight_model")
  save_model(tryCatch(fit_prior_wide_model(main_model), error = function(e) NULL),
             output_dir, "prior_wide_model")

  log_message("Finished: ", intervention)

}

for(intervention in selected_interventions()){
  tryCatch(run_sensitivity_models(intervention),
           error = function(e) log_message(intervention, ": ", e$message))
}

log_message("Sensitivity models complete.")
