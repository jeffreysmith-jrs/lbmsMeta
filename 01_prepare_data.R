#################################################
## 01_prepare_data.R
##
## Reads the database (data/observations.csv) and
## adds the moderator variables used in the models:
##   - sample sizes (N_int, N_ref)
##   - intervention age (intervention_length, years)
##   - intervention area (mean_size_ha, hectares)
##   - broad taxonomic group (taxa_grouped)
##   - biome (forested: Forested / Open) and climate
##     zone (zone: Tropical / Temperate)
##   - high_n_study: the top 10% of studies by number
##     of observations (sensitivity analysis)
##   - sampling effort (effort_int, effort_ref) and
##     summed_uneven_effort: summed richness with
##     unequal or unknown effort, left out of the models
##
## Reforestation and afforestation are recorded
## together as forestation. They are written as
## forestation.csv (pooled, for the forestation
## moderator models) and split by the author-defined
## forestation_type into reforestation.csv and
## aforestation.csv.
##
## Output: output/analysis_data/<intervention>.csv
#################################################

library(dplyr)
library(terra)

paths_file <- Find(file.exists, c("R/paths.R", file.path(Sys.getenv("LBMS_ANALYSIS_DIR"), "R", "paths.R")))
if(is.null(paths_file)) stop("Run from the repository root (the folder that contains run_analysis.R).")
source(paths_file)
source(file.path(code_dir, "R", "helpers.R"))

dir.create(analysis_data_dir, recursive = TRUE, showWarnings = FALSE)

#################################################
## Reference data
#################################################

## Surveyed taxon -> broad group
taxa_lookup <- read.csv(taxa_groups_file, stringsAsFactors = FALSE)

taxa_dict <- setNames(taxa_lookup$taxonomic_group, trimws(taxa_lookup$surveyed_taxon))

bad_groups <- setdiff(unique(taxa_lookup$taxonomic_group), taxa_groups)

if(length(bad_groups) > 0){
  stop("taxa_groups.csv has group names other than ", paste(taxa_groups, collapse = "/"), ": ",
       paste(bad_groups, collapse = ", "))
}

if(!file.exists(ecoregions_file)){
  stop("RESOLVE Ecoregions 2017 shapefile not found at ", ecoregions_file, " (see data/README.md).")
}

## RESOLVE Ecoregions 2017
ecoregions <- vect(ecoregions_file)

#################################################
## Ranges, age and area
#################################################

## Midpoint of a range written as text: "10-20" -> 15,
## "1995-2000" -> 1997.5, "30" -> 30. Anything that
## cannot be parsed is NA.
midpoint_from_range <- function(x){

  if(is.na(x)) return(NA_real_)

  x <- trimws(as.character(x))

  if(x == "") return(NA_real_)

  if(!grepl("-", x)) suppressWarnings(return(as.numeric(x)))

  vals <- suppressWarnings(as.numeric(trimws(strsplit(x, "-")[[1]])))
  vals <- vals[!is.na(vals)]

  if(length(vals) == 0) return(NA_real_)

  mean(vals)

}

midpoint_vector <- function(x) sapply(x, midpoint_from_range)

## Years between the (midpoint) year the intervention
## was established and the (midpoint) survey year.
## Negative ages can only come from data errors and
## are set to NA.
calculate_intervention_age <- function(intervention_year, survey_year){

  age <- midpoint_vector(survey_year) - midpoint_vector(intervention_year)
  age[!is.na(age) & age < 0] <- NA_real_
  age

}

## Area in hectares; NA (with a warning) for units
## that are not recognised.
convert_area <- function(area_value, area_units){

  if(is.na(area_value) || is.na(area_units)) return(NA_real_)

  area_value <- suppressWarnings(as.numeric(area_value))

  if(is.na(area_value)) return(NA_real_)

  area_units <- trimws(tolower(as.character(area_units)))

  if(area_units %in% c("hectare", "hectares", "ha",
                       "hectometers squared", "hectometres squared", "hm2", "hm^2")){
    return(area_value)
  }

  if(area_units %in% c("square kilometer", "square kilometers", "square kilometres",
                       "kilometers squared", "kilometres squared", "km2", "km^2",
                       "kilmeters squares", "kilometers squares")){
    return(area_value * 100)
  }

  if(area_units %in% c("square meter", "square meters", "square metre", "square metres",
                       "meters squared", "metres squared", "m2", "m^2",
                       "meters squares", "meters sqaured")){
    return(area_value / 10000)
  }

  if(area_units %in% c("acre", "acres")){
    return(area_value * 0.404686)
  }

  ## Area of a circle with this radius
  if(area_units %in% c("kilometer radius", "kilometre radius", "km radius")){
    return(pi * area_value^2 * 100)
  }

  if(area_units %in% c("centimeters squared", "centimetres squared", "cm2", "cm^2")){
    return(area_value / 1e8)
  }

  warning(paste("Unknown area unit:", area_units))

  NA_real_

}

calculate_area_hectares <- function(area_values, area_units){
  mapply(convert_area, midpoint_vector(area_values), area_units)
}

#################################################
## Moderator variables
#################################################

## Sample size per habitat: the larger of replicates
## per patch and number of patches.
add_sample_sizes <- function(df){

  df$N_int <- pmax(as.numeric(df$reps_per_patch_intervention),
                   as.numeric(df$patches_intervention), na.rm = TRUE)

  df$N_ref <- pmax(as.numeric(df$reps_per_patch_reference),
                   as.numeric(df$patches_reference), na.rm = TRUE)

  df

}

add_time_variables <- function(df){

  df$intervention_length <-
    calculate_intervention_age(df$year_intervention_applied, df$year_survey_start)

  df

}

add_area_variables <- function(df){

  df$mean_size_ha <- calculate_area_hectares(df$area_intervention, df$area_units)

  df

}

add_taxa_groups <- function(df){

  df$taxa_grouped <- taxa_dict[trimws(df$surveyed_taxon)]

  unmatched <- unique(df$surveyed_taxon[!is.na(df$surveyed_taxon) & is.na(df$taxa_grouped)])

  if(length(unmatched) > 0){
    log_message(
      "WARNING: surveyed_taxon values not in taxa_groups.csv (taxa_grouped = NA): ",
      paste(unmatched, collapse = "; ")
    )
  }

  df

}

## Biome from the site coordinates (RESOLVE Ecoregions
## 2017), and climate zone from absolute latitude.
##
## An observation that pools several sites carries their
## coordinates as a "; "-separated list. It is placed at
## the mean of its sites (map point and climate zone),
## and takes the biome (and ecoregion) that most of its
## sites fall in; sites off land are ignored, and a tie
## goes to the biome of the earliest-listed site.
##
## A site that falls in no ecoregion polygon because the
## polygons' coastline is coarse (coastal wetlands,
## estuaries, small islands, lake shores) takes the
## nearest ecoregion within coast_snap_km. An observation
## with no biome has forested = NA, so it is left out of
## the biome model rather than counted as open.
split_coords <- function(x){
  lapply(strsplit(as.character(x), ";", fixed = TRUE), safe_numeric)
}

snap_to_nearest_ecoregion <- function(site_eco, sites){

  attr_cols <- setdiff(names(site_eco), names(site_eco)[1])
  window <- coast_snap_km / 111 * 2

  for(i in which(is.na(site_eco$BIOME_NUM))){
    lon <- sites$lon[i]
    lat <- sites$lat[i]
    near <- tryCatch({
      box <- ext(lon - window / cos(lat * pi / 180), lon + window / cos(lat * pi / 180),
                 lat - window, lat + window)
      cand <- sf::st_make_valid(sf::st_as_sf(crop(ecoregions, box)))
      cand <- cand[!is.na(cand$BIOME_NUM), ]
      if(nrow(cand) == 0) NULL else {
        pt <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
        d_km <- as.numeric(sf::st_distance(pt, cand)) / 1000
        if(min(d_km) <= coast_snap_km) sf::st_drop_geometry(cand)[which.min(d_km), , drop = FALSE] else NULL
      }
    }, error = function(e) NULL)
    if(!is.null(near)) site_eco[i, attr_cols] <- near[1, attr_cols]
  }

  site_eco

}

add_geography <- function(df){

  lat_list <- split_coords(df$site_latitude_dd)
  lon_list <- split_coords(df$site_longitude_dd)
  n_sites <- pmax(lengths(lat_list), lengths(lon_list))
  usable <- lengths(lat_list) == lengths(lon_list) &
    vapply(lat_list, function(v) !anyNA(v), TRUE) &
    vapply(lon_list, function(v) !anyNA(v), TRUE)

  sites <- data.frame(
    row = rep(seq_len(nrow(df)), ifelse(usable, n_sites, 0)),
    lat = unlist(lat_list[usable]),
    lon = unlist(lon_list[usable])
  )

  pts <- vect(
    as.matrix(sites[, c("lon", "lat")]),
    type = "points",
    crs = "+proj=longlat +datum=WGS84"
  )
  site_eco <- terra::extract(ecoregions, pts)
  site_eco <- snap_to_nearest_ecoregion(site_eco, sites)

  ## One ecoregion row per observation: its single site,
  ## or the first site in the most common biome
  pick <- rep(NA_integer_, nrow(df))
  for(r in unique(sites$row)){
    idx <- which(sites$row == r)
    b <- site_eco$BIOME_NUM[idx]
    if(all(is.na(b))){
      pick[r] <- idx[1]
      next
    }
    counts <- table(factor(b[!is.na(b)], levels = unique(b[!is.na(b)])))
    pick[r] <- idx[which(b == names(counts)[which.max(counts)])[1]]
  }

  df$site_latitude_dd <- NA_real_
  df$site_longitude_dd <- NA_real_
  df$site_latitude_dd[usable] <- vapply(lat_list[usable], mean, 0)
  df$site_longitude_dd[usable] <- vapply(lon_list[usable], mean, 0)
  df$n_sites_coordinates <- ifelse(usable, n_sites, NA_integer_)

  eco_rows <- site_eco[pick, , drop = FALSE]
  eco_rows[[1]] <- seq_len(nrow(df))
  rownames(eco_rows) <- NULL
  df <- cbind(df, eco_rows)

  df$forested <- ifelse(is.na(df$BIOME_NUM), NA_character_,
                        ifelse(df$BIOME_NUM %in% forested_biome_ids, "Forested", "Open"))

  df$zone <- ifelse(abs(df$site_latitude_dd) < tropical_cutoff, "Tropical", "Temperate")

  df

}

## Sampling effort on each side = patches x replicates
## per patch x subsamples per replicate (components not
## recorded count as 1; NA when none is recorded).
## Summed richness grows with effort, so a summed
## observation is comparable only at equal effort;
## summed_uneven_effort flags summed observations whose
## effort differs between intervention and reference or
## is unknown. Means are per sampling unit and are never
## flagged.
sampling_effort <- function(...){
  parts <- lapply(list(...), function(x) suppressWarnings(as.numeric(x)))
  m <- do.call(cbind, parts)
  out <- apply(m, 1, function(r) if(all(is.na(r))) NA_real_ else prod(r, na.rm = TRUE))
  as.numeric(out)
}

add_effort_flag <- function(df){

  col <- function(name) if(name %in% names(df)) df[[name]] else rep(NA, nrow(df))

  df$effort_int <- sampling_effort(col("patches_intervention"), col("reps_per_patch_intervention"),
                                   col("subsamp_per_rep_intervention"))
  df$effort_ref <- sampling_effort(col("patches_reference"), col("reps_per_patch_reference"),
                                   col("subsamp_per_rep_reference"))

  summed <- trimws(as.character(df$S_aggregation)) %in% c("sum", "summed")
  equal <- !is.na(df$effort_int) & !is.na(df$effort_ref) &
    abs(df$effort_int - df$effort_ref) < 1e-9

  df$summed_uneven_effort <- summed & !equal

  df

}

## Flags every observation from the top 10% of
## studies by number of observations. Ties at the
## cut-off are all flagged.
add_high_n_flag <- function(df){

  study_counts <- as.data.frame(table(df$study_id), stringsAsFactors = FALSE)
  names(study_counts) <- c("study_id", "n_obs")

  n_top <- max(1, ceiling(0.10 * nrow(study_counts)))

  threshold <- sort(study_counts$n_obs, decreasing = TRUE)[n_top]

  df$high_n_study <- df$study_id %in% study_counts$study_id[study_counts$n_obs >= threshold]

  df

}

#################################################
## Read the database
##
## observations.csv holds every extracted comparison
## for all interventions; each intervention keeps the
## row order of the database.
#################################################

observations <- read.csv(observations_file, stringsAsFactors = FALSE, check.names = FALSE,
                         na.strings = c("", "NA"), encoding = "UTF-8")

read_intervention <- function(key){
  df <- observations[observations$intervention == key, setdiff(names(observations), "intervention")]
  rownames(df) <- NULL
  df
}

add_moderators <- function(df){

  df <- add_sample_sizes(df)
  df <- add_time_variables(df)
  df <- add_area_variables(df)
  df <- add_taxa_groups(df)
  df <- add_geography(df)
  df <- add_effort_flag(df)
  add_high_n_flag(df)

}

## Ecoregion attributes kept from RESOLVE Ecoregions 2017
ecoregion_columns <- c("ECO_NAME", "BIOME_NUM", "BIOME_NAME", "REALM", "ECO_ID")

write_analysis_data <- function(df, intervention){
  drop <- setdiff(names(df)[(match("n_sites_coordinates", names(df)) + 1):(match("forested", names(df)) - 1)],
                  ecoregion_columns)
  write.csv(df[, setdiff(names(df), drop)], file.path(analysis_data_dir, paste0(intervention, ".csv")),
            row.names = FALSE)
}

#################################################
## Interventions with one search each
#################################################

for(intervention in setdiff(unique(observations$intervention), "forestation")){

  log_message("Processing ", intervention)

  write_analysis_data(add_moderators(read_intervention(intervention)), intervention)

}

#################################################
## Forestation: pooled, then split into
## reforestation and afforestation. The high-N flag
## is recomputed within each split.
#################################################

log_message("Processing forestation")

forestation <- add_moderators(read_intervention("forestation"))

write_analysis_data(forestation, "forestation")

afforestation <- forestation[tolower(forestation$forestation_type) %in% c("aforestation", "afforestation"), ]

write_analysis_data(add_high_n_flag(afforestation), "aforestation")

reforestation <- forestation[tolower(forestation$forestation_type) == "reforestation", ]

write_analysis_data(add_high_n_flag(reforestation), "reforestation")

log_message("Forestation observations: ", nrow(forestation),
            " (reforestation ", nrow(reforestation), ", afforestation ", nrow(afforestation), ")")

#################################################
## Missing values in the model variables
#################################################

model_columns <- c(
  "study_id", "S_reference", "S_intervention", "S_aggregation",
  "taxa_grouped", "zone", "forested",
  "reference_habitat_structure_axis", "reference_habitat_human_use_axis",
  "experiment_type", "mean_size_ha", "intervention_length"
)

coverage <- do.call(rbind, lapply(c(interventions, "forestation"), function(intervention){

  d <- read.csv(file.path(analysis_data_dir, paste0(intervention, ".csv")), stringsAsFactors = FALSE)

  na_counts <- vapply(model_columns, function(col){
    if(!col %in% names(d)) return(NA_integer_)
    as.integer(sum(is.na(d[[col]]) | trimws(as.character(d[[col]])) %in% c("", "NA")))
  }, integer(1))

  data.frame(intervention = intervention, rows = nrow(d),
             t(setNames(na_counts, paste0("na_", model_columns))), check.names = FALSE)

}))

write.csv(coverage, file.path(analysis_data_dir, "_coverage_model_columns.csv"), row.names = FALSE)

log_message("Missing values in model columns (rows / NAs per column):")
print(coverage, row.names = FALSE)
