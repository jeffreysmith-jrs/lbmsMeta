#################################################
## 03_fit_main_models.R
##
## Bayesian meta-analysis models (brms + cmdstanr)
## for each intervention:
##   - main model: lrr_richness ~ 1 + (1 | study_id)
##     + S_aggregation
##   - the main model plus one moderator: taxonomic
##     group, climate zone, biome, reference habitat
##     structure (woody/open), reference habitat human
##     use, study design, log(area), and age.
## For the pooled forestation data (reforestation +
## afforestation): biome, zone, reference habitat
## structure and human use, study design, and
## forestation method x biome.
##
## Rows whose S_aggregation is neither mean (medians
## count as means) nor aggregated, and summed rows with
## unequal sampling effort, are dropped
## (modelled_rows(), R/models.R), and S_aggregation is
## sum-to-zero coded (R/models.R). Each moderator
## model also drops rows missing that moderator.
##
## Output: output/models/<intervention>/<model>.rds
## (and main_draws.csv, the main model's posterior
## draws).
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

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

read_effect_data <- function(intervention){

  df <- read.csv(file.path(effects_dir, paste0(intervention, ".csv")), stringsAsFactors = FALSE)
  df <- df[modelled_rows(df), ]

  df$S_aggregation <- simplify_aggregation(df$S_aggregation)
  df$S_aggregation <- sum_code_aggregation(df$S_aggregation)

  df[!is.na(df$lrr_richness), ]

}

## Shared brm() arguments. The log rows are labelled
## with the fit_*_model() function that called this.
fit_brm <- function(formula, df, prior = priors_moderator){

  brm(
    formula,
    data = df,
    family = gaussian(),
    prior = prior,
    chains = n_chains,
    iter = n_iter,
    warmup = n_warmup,
    cores = n_cores,
    control = brms_control_default,
    backend = "cmdstanr",
    seed = project_seed,
    caller = deparse(sys.call(-1)[[1]])
  )

}

## Factor with the reference level first, using only
## the levels present in this intervention's data
relevel_present <- function(x, preferred_order){
  factor(x, levels = intersect(preferred_order, unique(x)))
}

#################################################
## Model functions
#################################################

fit_main_model <- function(df){
  fit_brm(formula_main, df)
}

## A moderator model needs at least two levels of
## its moderator; otherwise it is skipped (NULL).
has_two_levels <- function(x) length(unique(na.omit(as.character(x)))) > 1

fit_taxa_model <- function(df){
  if(!has_two_levels(df$taxa_grouped)) return(NULL)
  df$taxa_grouped <- factor(df$taxa_grouped)
  fit_brm(formula_taxa, df)
}

fit_zone_model <- function(df){
  if(!has_two_levels(df$zone)) return(NULL)
  fit_brm(formula_zone, df)
}

fit_biome_model <- function(df){
  if(!has_two_levels(df$forested)) return(NULL)
  fit_brm(formula_biome, df)
}

## Woody is the reference level
fit_structure_model <- function(df){
  df <- df[!is.na(df$reference_habitat_structure_axis), ]
  if(!has_two_levels(df$reference_habitat_structure_axis)) return(NULL)
  df$reference_habitat_structure_axis <-
    relevel_present(df$reference_habitat_structure_axis, c("Woody", "Open"))
  fit_brm(formula_structure, df)
}

## Intact natural is the reference level
fit_human_use_model <- function(df){
  df <- df[!is.na(df$reference_habitat_human_use_axis), ]
  if(!has_two_levels(df$reference_habitat_human_use_axis)) return(NULL)
  df$reference_habitat_human_use_axis <-
    relevel_present(df$reference_habitat_human_use_axis, c("Intact natural", "Degraded natural", "Human use"))
  fit_brm(formula_human_use, df)
}

## Natural (observational) is the reference level
fit_experiment_type_model <- function(df){
  df <- df[!is.na(df$experiment_type), ]
  df$experiment_type <- relevel_present(df$experiment_type, c("Natural", "Manipulative"))
  fit_brm(formula_experiment_type, df)
}

fit_area_model <- function(df){
  fit_brm(formula_area, df[!is.na(df$mean_size_ha), ])
}

fit_time_model <- function(df){
  fit_brm(formula_time, df[!is.na(df$intervention_length), ])
}

## Natural forest restoration vs timber plantation,
## in historically forested vs open biomes. The other
## methods (riparian corridors, windbreaks, fruit and
## nut plantations) have too few comparisons.
fit_forestation_type_model <- function(df){
  df <- df[df$forest_type %in% c("natural forest restoration", "timber plantation") &
             !is.na(df$forest_type) & !is.na(df$forested), ]
  df$forest_type <- factor(df$forest_type)
  df$forested <- factor(df$forested)
  fit_brm(formula_forestation_method, df)
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

## Study design needs both levels present (solar and
## wind are entirely observational).
has_two_experiment_types <- function(df){
  "experiment_type" %in% names(df) && length(unique(na.omit(df$experiment_type))) > 1
}

run_intervention_models <- function(intervention){

  log_message("Running ", intervention)

  df <- read_effect_data(intervention)
  set_model_base(intervention, df)

  output_dir <- file.path(results_dir, intervention)
  dir.create(output_dir, showWarnings = FALSE)

  m_main <- fit_main_model(df)
  save_model(m_main, output_dir, "main_model")
  if(!is.null(m_main)){
    write.csv(posterior::as_draws_df(m_main), file.path(output_dir, "main_draws.csv"), row.names = FALSE)
  }

  save_model(fit_taxa_model(df), output_dir, "taxa_model")
  save_model(fit_zone_model(df), output_dir, "zone_model")
  save_model(fit_biome_model(df), output_dir, "biome_model")
  save_model(fit_structure_model(df), output_dir, "structure_model")
  save_model(fit_human_use_model(df), output_dir, "human_use_model")

  if(has_two_experiment_types(df)){
    save_model(fit_experiment_type_model(df), output_dir, "experiment_type_model")
  } else {
    log_message(intervention, ": skipping experiment_type_model (one study design only)")
  }

  save_model(fit_area_model(df), output_dir, "area_model")
  save_model(fit_time_model(df), output_dir, "time_model")

}

run_forestation_models <- function(){

  log_message("Running forestation")

  df <- read_effect_data("forestation")
  set_model_base("forestation", df)

  output_dir <- file.path(results_dir, "forestation")
  dir.create(output_dir, showWarnings = FALSE)

  save_model(fit_biome_model(df), output_dir, "biome_model")
  save_model(fit_zone_model(df), output_dir, "zone_model")
  save_model(fit_structure_model(df), output_dir, "structure_model")
  save_model(fit_human_use_model(df), output_dir, "human_use_model")

  if(has_two_experiment_types(df)){
    save_model(fit_experiment_type_model(df), output_dir, "experiment_type_model")
  }

  save_model(fit_forestation_type_model(df), output_dir, "forestation_type_model")

}

for(intervention in selected_interventions()){
  tryCatch(run_intervention_models(intervention),
           error = function(e) log_message(intervention, ": ", e$message))
}

if("forestation" %in% selected_interventions(c(interventions, "forestation"))){
  tryCatch(run_forestation_models(),
           error = function(e) log_message("forestation: ", e$message))
}

log_message("Main and moderator models complete.")
