#################################################
## helpers.R
##
## Constants and small utilities shared by every
## step of the analysis.
#################################################

#################################################
## Interventions
##
## Internal keys (file and folder names) and the
## labels used in figures and tables. "aforestation"
## follows the spelling of the source data.
## "forestation" is reforestation + afforestation
## pooled, used for the biome, climate-zone,
## reference-habitat, study-design and forestation-
## method moderator models.
#################################################

interventions <- c(
  "reforestation",
  "grassland_restoration",
  "wetland_restoration-inland",
  "wetland_restoration-coastal",
  "agroforestry",
  "aforestation",
  "beccs",
  "solar",
  "wind"
)

plot_labels <- c(
  "reforestation" = "Reforestation",
  "grassland_restoration" = "Grassland restoration",
  "wetland_restoration-inland" = "Inland wetland restoration",
  "wetland_restoration-coastal" = "Coastal wetland restoration",
  "agroforestry" = "Agroforestry",
  "aforestation" = "Afforestation",
  "beccs" = "Bioenergy cropping",
  "solar" = "Solar",
  "wind" = "Wind",
  "forestation" = "Forestation"
)

#################################################
## Classification constants (Methods)
#################################################

## Broad taxonomic groups (taxa_groups.csv)
taxa_groups <- c("Plants", "Invertebrates", "Vertebrates", "Microbes")

## Tropical: absolute latitude < 23.5 degrees
tropical_cutoff <- 23.5

## RESOLVE Ecoregions 2017 BIOME_NUM values counted as
## historically forested; every other biome is open.
forested_biome_ids <- c(
  1,   # Tropical & subtropical moist broadleaf forests
  2,   # Tropical & subtropical dry broadleaf forests
  3,   # Tropical & subtropical coniferous forests
  4,   # Temperate broadleaf & mixed forests
  5,   # Temperate conifer forests
  6,   # Boreal forests/taiga
  14   # Mangroves
)

## A site that falls in no ecoregion polygon (coastal
## marsh, estuary, small island, lake shore) takes the
## nearest ecoregion if it is within this distance (km)
coast_snap_km <- 15

## Seed for every model fit and every random draw
project_seed <- 34953

#################################################
## Utilities
#################################################

log_message <- function(...){
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
}

## as.numeric() that turns anything unparseable into NA
safe_numeric <- function(x){
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

## The interventions a model step fits: all of `all`,
## or those named in the comma-separated environment
## variable LBMS_INTERVENTIONS (e.g.
## LBMS_INTERVENTIONS=beccs,solar), to refit a subset
selected_interventions <- function(all = interventions){
  chosen <- trimws(strsplit(Sys.getenv("LBMS_INTERVENTIONS"), ",")[[1]])
  if(length(chosen) == 0) return(all)
  unknown <- setdiff(chosen, c(interventions, "forestation"))
  if(length(unknown) > 0) stop("Unknown intervention(s) in LBMS_INTERVENTIONS: ", paste(unknown, collapse = ", "))
  intersect(all, chosen)
}
