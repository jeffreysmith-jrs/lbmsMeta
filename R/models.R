#################################################
## models.R
##
## Model settings shared by steps 03-06: sampler
## settings, priors, model formulas, the coding of
## S_aggregation, and the wrapper every model is
## fitted through (row-count log + convergence
## check with the refit sequence in the Methods).
#################################################

library(brms)

set.seed(project_seed)

#################################################
## Sampler settings
##
## Default: 4 chains, 5,000 iterations (2,000
## warm-up), adapt_delta 0.99, max tree depth 15.
## Long run, used when a model has not converged:
## 7,500 iterations (2,500 warm-up), adapt_delta
## 0.999, max tree depth 20. Extra-long run: 15,000
## iterations (5,000 warm-up), same sampler settings.
## Final refit, for divergent transitions that
## remain: the long run with adapt_delta 0.99999
## (smaller step size).
#################################################

n_chains <- 4
n_cores <- 4
n_iter <- 5000
n_warmup <- 2000
brms_control_default <- list(adapt_delta = 0.99, max_treedepth = 15)

n_iter_long <- 7500
n_warmup_long <- 2500
brms_control_long <- list(adapt_delta = 0.999, max_treedepth = 20)

n_iter_extra_long <- 15000
n_warmup_extra_long <- 5000

## Final refit for a model that still has not converged
brms_control_final <- list(adapt_delta = 0.99999, max_treedepth = 20)

#################################################
## Priors
##
## priors_moderator: every model with fixed effects,
##   including the main model (whose S_aggregation
##   term is a fixed effect).
## priors_intercept: intercept-only models with a
##   study random intercept (sensitivity analyses).
## priors_tight / priors_wide: the prior-sensitivity
##   refits of the main model. The three schemes
##   differ only in scale (Intercept and b: 1 / 5 /
##   10; between-study SD: 0.5 / 1 / 2).
## The residual SD (sigma) keeps the brms default,
## half-Student-t(3, 0, 2.5), throughout.
#################################################

priors_intercept <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(cauchy(0, 1), class = "sd")
)

priors_moderator <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 5), class = "b"),
  prior(cauchy(0, 1), class = "sd")
)

priors_tight <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 0.5), class = "sd")
)

priors_wide <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 10), class = "b"),
  prior(cauchy(0, 2), class = "sd")
)

#################################################
## Model formulas
##
## Every model regresses the log response ratio of
## species richness on its moderator(s), with a
## random intercept for study and a fixed effect of
## aggregation method (mean vs aggregated richness).
#################################################

formula_main <- bf(lrr_richness ~ 1 + (1 | study_id) + S_aggregation)

formula_taxa <- bf(lrr_richness ~ taxa_grouped + (1 | study_id) + S_aggregation)

## Tropical vs temperate
formula_zone <- bf(lrr_richness ~ zone + (1 | study_id) + S_aggregation)

## Historically forested vs open biome
formula_biome <- bf(lrr_richness ~ forested + (1 | study_id) + S_aggregation)

## Reference habitat land cover: Woody vs Open
formula_structure <- bf(lrr_richness ~ reference_habitat_structure_axis + (1 | study_id) + S_aggregation)

## Reference habitat land use: Intact natural /
## Degraded natural / Human use
formula_human_use <- bf(lrr_richness ~ reference_habitat_human_use_axis + (1 | study_id) + S_aggregation)

## Study design: Natural (observational) vs Manipulative
formula_experiment_type <- bf(lrr_richness ~ experiment_type + (1 | study_id) + S_aggregation)

## Intervention area, log(hectares)
formula_area <- bf(lrr_richness ~ log(mean_size_ha) + (1 | study_id) + S_aggregation)

## Years since the intervention was established
formula_time <- bf(lrr_richness ~ intervention_length + (1 | study_id) + S_aggregation)

## Forestation method x biome (pooled forestation data)
formula_forestation_method <- bf(lrr_richness ~ forest_type * forested + S_aggregation + (1 | study_id))

#################################################
## Aggregation method (S_aggregation)
##
## Extraction recorded several spellings. They are
## collapsed to "mean" (mean richness across plots,
## with or without an uncertainty term) and
## "aggregated" (richness summed across plots).
## Medians are treated as means (for box plots, step
## 03 takes the median as the mean and the
## interquartile range as the uncertainty), so they
## are also "mean". Other values (Chao, rarefied,
## ...) are returned as NA and those rows are dropped
## from the models.
#################################################

simplify_aggregation <- function(x){

  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("mean", "mean ± uncert", "median", "median ± uncert")] <- "mean"
  out[x %in% c("sum", "summed")] <- "aggregated"
  out

}

## The rows every model (and every count of modelled
## observations) uses: aggregation method mean or
## aggregated (simplify_aggregation()), and not summed
## richness with unequal or unknown sampling effort
## between intervention and reference
## (summed_uneven_effort, step 01). Summed richness
## grows with sampling effort, so those comparisons are
## not valid; like Chao and rarefied estimates, the
## rows stay in the data files but no model uses them.
## Rows without an effect size are dropped by each
## model separately.
modelled_rows <- function(df){
  !is.na(simplify_aggregation(df$S_aggregation)) & !(df$summed_uneven_effort %in% TRUE)
}

## Sum-to-zero contrasts, so the model intercept is
## the average over the two aggregation methods rather
## than the value for one of them.
sum_code_aggregation <- function(x){

  x <- factor(x)
  contrasts(x) <- contr.sum(nlevels(x))
  x

}

#################################################
## Model fitting wrapper
##
## Steps 04-05 call brm() through brm_logged(), which
##   1. appends to output/models/model_na_drops.csv how
##      many rows the model uses: rows in the
##      intervention's data (set_model_base()), rows
##      passed to brm() after the model's own filters,
##      and complete cases on the formula variables
##      (what brms fits);
##   2. checks convergence and, if needed, refits
##      following the Methods (ensure_converged()).
## Each log row is labelled with `caller`, by default
## the name of the function that called brm().
#################################################

model_log_env <- new.env()


set_model_base <- function(intervention, df){

  model_log_env$intervention <- intervention
  model_log_env$n_base <- nrow(df)
  invisible(NULL)

}

append_log_row <- function(row, file_name){

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(results_dir, file_name)
  write.table(row, out, sep = ",", row.names = FALSE,
              col.names = !file.exists(out), append = file.exists(out))

}

log_model_na_drops <- function(formula, data, caller){

  f <- if(inherits(formula, "brmsformula")) formula$formula else formula
  vars <- intersect(all.vars(f), names(data))
  n_used <- sum(stats::complete.cases(data[, vars, drop = FALSE]))
  n_base <- if(is.null(model_log_env$n_base)) NA_integer_ else model_log_env$n_base

  row <- data.frame(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    intervention = if(is.null(model_log_env$intervention)) NA_character_ else model_log_env$intervention,
    model = caller,
    n_intervention = n_base,
    n_passed_to_brm = nrow(data),
    n_used = n_used,
    dropped_before_brm = n_base - nrow(data),
    dropped_na_in_brm = nrow(data) - n_used,
    variables = paste(vars, collapse = " + "),
    stringsAsFactors = FALSE
  )

  log_message(
    row$intervention, " / ", caller, ": ", n_used, " of ", n_base,
    " rows used (", row$dropped_before_brm, " dropped before brm, ",
    row$dropped_na_in_brm, " dropped by brm for NA in: ", row$variables, ")"
  )

  append_log_row(row, "model_na_drops.csv")
  invisible(row)

}

brm_logged <- function(formula, data, ..., caller = NULL){

  if(is.null(caller)){
    caller <- tryCatch(deparse(sys.call(-1)[[1]]), error = function(e) NA_character_)
  }

  log_model_na_drops(formula, data, caller)

  fit <- brms::brm(formula = formula, data = data, ...)

  ensure_converged(fit, caller)

}

#################################################
## Convergence check and refit sequence
##
## A model has converged when max R-hat <= 1.01,
## min bulk and tail ESS >= 400 and there are no
## divergent transitions. If not:
##   1. refit with the long-run sampler settings;
##   2. if still not converged, refit with the
##      long-run settings, seed + 1 and every
##      parameter initialised at 0;
##   3. if still not converged, refit with the
##      extra-long run (original seed);
##   4. if still not converged, refit with the
##      long-run iterations and adapt_delta 0.99999;
##   5. if still not converged, keep that last fit and
##      flag it (it is reported as not converged in
##      the convergence table).
## Formula, data and priors never change. Every
## attempt is appended to
## output/models/model_convergence_log.csv.
## Leave-one-out refits (step 04) are not checked.
#################################################

rhat_max <- 1.01

ess_min <- 400

model_convergence_diagnostics <- function(fit){

  ess <- posterior::summarise_draws(posterior::as_draws_array(fit), "ess_bulk", "ess_tail")
  np <- brms::nuts_params(fit)

  d <- list(
    max_rhat = max(brms::rhat(fit), na.rm = TRUE),
    min_ess_bulk = min(ess$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(ess$ess_tail, na.rm = TRUE),
    n_divergent = sum(np$Value[np$Parameter == "divergent__"])
  )

  d$converged <-
    d$max_rhat <= rhat_max &&
    d$min_ess_bulk >= ess_min &&
    d$min_ess_tail >= ess_min &&
    d$n_divergent == 0

  d

}

log_model_convergence <- function(caller, attempt, d, status){

  row <- data.frame(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    intervention = if(is.null(model_log_env$intervention)) NA_character_ else model_log_env$intervention,
    model = caller,
    attempt = attempt,
    max_rhat = if(is.null(d)) NA else round(d$max_rhat, 4),
    min_ess_bulk = if(is.null(d)) NA else round(d$min_ess_bulk),
    min_ess_tail = if(is.null(d)) NA else round(d$min_ess_tail),
    n_divergent = if(is.null(d)) NA else d$n_divergent,
    status = status,
    stringsAsFactors = FALSE
  )

  log_message(
    row$intervention, " / ", caller, " [", attempt, "]: ", status,
    " (max R-hat ", row$max_rhat, ", min ESS bulk/tail ", row$min_ess_bulk,
    "/", row$min_ess_tail, ", divergent ", row$n_divergent, ")"
  )

  append_log_row(row, "model_convergence_log.csv")
  invisible(row)

}

ensure_converged <- function(fit, caller = NA_character_){

  d <- model_convergence_diagnostics(fit)

  if(d$converged){
    log_model_convergence(caller, "0 default", d, "converged")
    return(fit)
  }

  log_model_convergence(caller, "0 default", d, "not converged")

  ## recompile = TRUE: a fit from another R session has
  ## no compiled Stan program attached.
  refit_steps <- list(
    "1 long-run" = list(seed = project_seed, init = NULL, iter = n_iter_long, warmup = n_warmup_long),
    "2 long-run, new seed, init 0" = list(seed = project_seed + 1, init = 0, iter = n_iter_long, warmup = n_warmup_long),
    "3 extra-long run" = list(seed = project_seed, init = NULL, iter = n_iter_extra_long, warmup = n_warmup_extra_long),
    "4 long-run, adapt_delta 0.99999" = list(seed = project_seed, init = NULL, iter = n_iter_long, warmup = n_warmup_long,
                                             control = brms_control_final)
  )

  last_fit <- fit

  for(attempt in names(refit_steps)){

    step <- refit_steps[[attempt]]

    args <- list(
      fit,
      recompile = TRUE,
      iter = step$iter,
      warmup = step$warmup,
      control = if(is.null(step$control)) brms_control_long else step$control,
      chains = n_chains,
      cores = n_cores,
      seed = step$seed
    )
    if(!is.null(step$init)) args$init <- step$init

    refit <- tryCatch(
      do.call(update, args),
      error = function(e){
        log_message("Refit (", attempt, ") failed for ", caller, ": ", conditionMessage(e))
        NULL
      }
    )

    if(is.null(refit)){
      log_model_convergence(caller, attempt, NULL, "refit error")
      next
    }

    d <- model_convergence_diagnostics(refit)
    last_fit <- refit

    if(d$converged){
      log_model_convergence(caller, attempt, d, "converged")
      return(refit)
    }

    log_model_convergence(caller, attempt, d, "not converged")

  }

  log_model_convergence(caller, "5 kept", model_convergence_diagnostics(last_fit),
                        "KEPT - not converged after all refits (flagged)")
  warning(caller, ": not converged after all refits - kept and flagged (see model_convergence_log.csv).")

  last_fit

}
