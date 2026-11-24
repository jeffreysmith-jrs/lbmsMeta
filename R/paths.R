#################################################
## paths.R
##
## Every folder and file the analysis reads or
## writes. Sourced at the top of every script.
##
## The repository root is the folder that contains
## run_analysis.R. It is found from the location of
## this file, from the working directory, or from the
## LBMS_ANALYSIS_DIR environment variable if set.
##
## Inputs are read from data/ (the Zenodo deposit,
## see data/README.md). Everything the analysis
## writes goes to output/, or to the folder named by
## LBMS_OUTPUT_DIR.
#################################################

find_repository_dir <- function(){

  env_dir <- Sys.getenv("LBMS_ANALYSIS_DIR")
  if(nzchar(env_dir)) return(normalizePath(env_dir, mustWork = TRUE))

  ## source() records the file it is reading as `ofile`
  for(frame in rev(sys.frames())){
    ofile <- frame$ofile
    if(!is.null(ofile) && basename(ofile) == "paths.R" && file.exists(ofile)){
      return(normalizePath(file.path(dirname(ofile), ".."), mustWork = TRUE))
    }
  }

  d <- normalizePath(getwd(), mustWork = TRUE)
  repeat{
    if(file.exists(file.path(d, "run_analysis.R"))) return(d)
    if(dirname(d) == d){
      stop("Could not find the repository root (the folder that contains run_analysis.R). ",
           "Run from inside it or set LBMS_ANALYSIS_DIR.")
    }
    d <- dirname(d)
  }

}

code_dir <- find_repository_dir()
rm(find_repository_dir)

output_dir <-
  if(nzchar(Sys.getenv("LBMS_OUTPUT_DIR"))){
    normalizePath(Sys.getenv("LBMS_OUTPUT_DIR"), mustWork = FALSE)
  } else {
    file.path(code_dir, "output")
  }

#################################################
## Inputs (data/)
#################################################

data_dir <- file.path(code_dir, "data")

## One row per extracted pair-wise comparison
observations_file <- file.path(data_dir, "observations.csv")

## Every study assessed at full text, with its outcome
studies_file <- file.path(data_dir, "studies.csv")

## Surveyed taxon -> broad taxonomic group
taxa_groups_file <- file.path(data_dir, "taxa_groups.csv")

## Full-text language of every paper assessed at full
## text, and Web of Science language of every record
## screened at title/abstract
paper_language_file <- file.path(data_dir, "paper_language.csv")
title_abstract_language_file <- file.path(data_dir, "title_abstract_language.csv")

## RESOLVE Ecoregions 2017 (Dinerstein et al. 2017)
ecoregions_file <- file.path(data_dir, "Ecoregions2017", "Ecoregions2017.shp")

#################################################
## Outputs (output/)
#################################################

analysis_data_dir <- file.path(output_dir, "analysis_data")  # 01
effects_dir       <- file.path(output_dir, "effects")        # 02
results_dir       <- file.path(output_dir, "models")         # 03-04 (fitted models), logs of every fit
summary_dir       <- file.path(output_dir, "summaries")      # 05, 06
figure_dir        <- file.path(output_dir, "figures")        # 06
table_dir         <- file.path(output_dir, "tables")         # 06
log_dir           <- file.path(output_dir, "logs")           # run_analysis.R
