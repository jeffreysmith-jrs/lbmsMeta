#################################################
## 02_effect_sizes.R
##
## Effect size for every observation: the log
## response ratio of species richness,
##   lrr_richness = log(S_intervention / S_reference),
## and, where the study reports mean richness with an
## uncertainty term, its standard error
## (lrr_se_richness) and variance (lrr_var_richness),
## used by the variance-weighted and publication-bias
## sensitivity analyses.
##
## Input:  output/analysis_data/<intervention>.csv
## Output: effects/<intervention>.csv
#################################################

paths_file <- Find(file.exists, c("R/paths.R", file.path(Sys.getenv("LBMS_ANALYSIS_DIR"), "R", "paths.R")))
if(is.null(paths_file)) stop("Run from the repository root (the folder that contains run_analysis.R).")
source(paths_file)
source(file.path(code_dir, "R", "helpers.R"))

dir.create(effects_dir, recursive = TRUE, showWarnings = FALSE)

numeric_cols <- c(
  "S_reference", "S_intervention",
  "N_ref", "N_int",
  "S_uncert_reference", "S_uncert_intervention",
  "Min_reference", "Lower_quartile_reference", "Median_reference",
  "Upper_quartile_reference", "Max_reference",
  "Min_intervention", "Lower_quartile_intervention", "Median_intervention",
  "Upper_quartile_intervention", "Max_intervention"
)

#################################################
## Reported uncertainty -> sample variance
##
## SD: squared. SE: multiplied by sqrt(n), then
## squared. 95% CI half-width: divided by 1.96,
## multiplied by sqrt(n), then squared. Box plot:
## the interquartile range divided by 1.35, then
## squared (Luo et al. 2018). Other types (Chao,
## "other") have no usable variance.
#################################################

uncert_to_s2 <- function(uncert, type, n = NA){

  uncert <- suppressWarnings(as.numeric(uncert))

  if(is.na(type)) return(NA_real_)

  type <- trimws(as.character(type))

  ## Spelling variants used during extraction
  type <- switch(type,
    "boxxplot" = "boxplot",
    "SE" = "SE (standard error)",
    "NS" = "other (explain in notes)",
    "SD" = "SD (σ/s/standard deviation)",
    type
  )

  if(type == "SD (σ/s/standard deviation)") return(uncert^2)

  if(type == "SE (standard error)"){
    if(is.na(n)) return(NA_real_)
    return((uncert * sqrt(n))^2)
  }

  if(type == "95% CI"){
    if(is.na(n)) return(NA_real_)
    return(((uncert / 1.96) * sqrt(n))^2)
  }

  if(type == "boxplot") return((uncert / 1.35)^2)

  if(type %in% c("chao", "other (explain in notes)")) return(NA_real_)

  warning(paste("Unknown uncertainty type:", type))

  NA_real_

}

#################################################
## Box-plot observations
##
## Observations extracted from box plots carry
## min / Q1 / median / Q3 / max instead of a mean.
## Where S_reference or S_intervention is missing,
## the mean is estimated as the median (or the
## midpoint of Q1 and Q3 if no median was recorded),
## and the uncertainty as the interquartile range
## Q3 - Q1 (converted to a variance in
## uncert_to_s2()). Values already entered are never
## overwritten.
#################################################

boxplot_mean <- function(median = NA, q1 = NA, q3 = NA){

  median <- suppressWarnings(as.numeric(median))
  q1 <- suppressWarnings(as.numeric(q1))
  q3 <- suppressWarnings(as.numeric(q3))

  if(!is.na(median)) median
  else if(!is.na(q1) & !is.na(q3)) (q1 + q3) / 2
  else NA_real_

}

impute_boxplot_summary_stats <- function(df){

  is_boxplot <-
    !is.na(df$uncert_type) &
    trimws(as.character(df$uncert_type)) %in% c("boxplot", "boxxplot")

  for(side in c("reference", "intervention")){

    mean_col <- paste0("S_", side)
    uncert_col <- paste0("S_uncert_", side)
    q1 <- df[[paste0("Lower_quartile_", side)]]
    q3 <- df[[paste0("Upper_quartile_", side)]]

    fill_mean <- is_boxplot & is.na(df[[mean_col]])
    if(any(fill_mean)){
      df[[mean_col]][fill_mean] <- mapply(
        boxplot_mean,
        median = df[[paste0("Median_", side)]][fill_mean],
        q1 = q1[fill_mean],
        q3 = q3[fill_mean]
      )
    }

    fill_uncert <- is_boxplot & is.na(df[[uncert_col]])
    if(any(fill_uncert)){
      df[[uncert_col]][fill_uncert] <- q3[fill_uncert] - q1[fill_uncert]
    }

  }

  df

}

#################################################
## Log response ratio and its standard error
## (Hedges et al. 1999)
#################################################

lrr <- function(xbar_int, xbar_ref){

  xbar_int <- as.numeric(xbar_int)
  xbar_ref <- as.numeric(xbar_ref)

  if(any(is.na(c(xbar_int, xbar_ref)))) return(NA_real_)
  if(xbar_int <= 0 || xbar_ref <= 0) return(NA_real_)

  log(xbar_int / xbar_ref)

}

lrr_se <- function(xbar_int, xbar_ref, n_int, n_ref, var_int, var_ref){

  if(any(is.na(c(xbar_int, xbar_ref, n_int, n_ref, var_int, var_ref)))) return(NA_real_)

  sqrt((var_int / (n_int * xbar_int^2)) + (var_ref / (n_ref * xbar_ref^2)))

}

#################################################
## Run
#################################################

process_effect_sizes <- function(intervention){

  log_message("Processing ", intervention)

  df <- read.csv(file.path(analysis_data_dir, paste0(intervention, ".csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)

  for(col in intersect(numeric_cols, names(df))){
    df[[col]] <- safe_numeric(df[[col]])
  }

  df <- impute_boxplot_summary_stats(df)

  df$var_ref <- mapply(uncert_to_s2, uncert = df$S_uncert_reference, type = df$uncert_type, n = df$N_ref)
  df$var_int <- mapply(uncert_to_s2, uncert = df$S_uncert_intervention, type = df$uncert_type, n = df$N_int)

  df$lrr_richness <- mapply(lrr, df$S_intervention, df$S_reference)

  df$lrr_se_richness <- mapply(
    lrr_se,
    xbar_int = df$S_intervention, xbar_ref = df$S_reference,
    n_int = df$N_int, n_ref = df$N_ref,
    var_int = df$var_int, var_ref = df$var_ref
  )

  df$lrr_var_richness <- df$lrr_se_richness^2

  log_message(intervention, ": ", nrow(df), " rows; ",
              sum(is.na(df$lrr_richness)), " without lrr_richness; ",
              sum(is.na(df$lrr_var_richness)), " without a variance")

  write.csv(df, file.path(effects_dir, paste0(intervention, ".csv")), row.names = FALSE)

}

for(intervention in c(interventions, "forestation")){
  process_effect_sizes(intervention)
}

log_message("Effect sizes complete.")
