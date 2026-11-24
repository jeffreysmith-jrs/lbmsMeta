#################################################
## run_analysis.R
##
## Runs the analysis steps in order, each in a fresh
## R process, and stops at the first failure. Each
## step's output is also written to
## output/logs/<step>_<timestamp>.log.
##
## From the repository root:
##   Rscript run_analysis.R                  # all steps
##   Rscript run_analysis.R --from 03        # 03 to 06
##   Rscript run_analysis.R --only 01,02,06  # figures and tables from the deposited model summaries
##   Rscript run_analysis.R --skip 04        # all but 04
##   Rscript run_analysis.R --dry-run        # list the steps only
##
## Approximate run times on a laptop: 01-02 about a
## minute; 03 and 04 several hours each; 05 about
## 5 minutes; 06 about a minute.
#################################################

paths_file <- Find(file.exists, c("R/paths.R", file.path(Sys.getenv("LBMS_ANALYSIS_DIR"), "R", "paths.R")))
if(is.null(paths_file)) stop("Run from the repository root (the folder that contains run_analysis.R).")
source(paths_file)

step_scripts <- c(
  "01" = "01_prepare_data.R",
  "02" = "02_effect_sizes.R",
  "03" = "03_fit_main_models.R",
  "04" = "04_fit_sensitivity_models.R",
  "05" = "05_summarise_models.R",
  "06" = "06_figures_and_tables.R"
)

#################################################
## Command-line options
#################################################

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag){
  i <- match(flag, args)
  if(is.na(i) || i == length(args)) return(NULL)
  steps <- trimws(strsplit(args[i + 1], ",")[[1]])
  bad <- setdiff(steps, names(step_scripts))
  if(length(bad) > 0){
    stop("Unknown step(s): ", paste(bad, collapse = ", "), ". Valid: ", paste(names(step_scripts), collapse = " "))
  }
  steps
}

run_steps <- setNames(rep(TRUE, length(step_scripts)), names(step_scripts))

if(!is.null(only <- arg_value("--only"))) run_steps[] <- names(run_steps) %in% only
if(!is.null(from <- arg_value("--from"))) run_steps[] <- seq_along(run_steps) >= match(from[1], names(run_steps))
if(!is.null(skip <- arg_value("--skip"))) run_steps[skip] <- FALSE

todo <- names(run_steps)[run_steps]
message("Steps to run: ", if(length(todo)) paste(todo, collapse = " ") else "(none)")
if("--dry-run" %in% args) quit(save = "no")

#################################################
## Run
#################################################

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
timings <- data.frame(step = character(0), minutes = numeric(0), status = character(0))

for(step in todo){

  script <- file.path(code_dir, step_scripts[[step]])
  log_file <- file.path(log_dir, paste0(step, "_", stamp, ".log"))
  message("\n==== ", step, ": ", basename(script), "  (log: ", log_file, ")")

  t0 <- Sys.time()
  cmd <- sprintf("cd %s && %s %s 2>&1 | tee %s; exit ${PIPESTATUS[0]}",
                 shQuote(code_dir), shQuote(rscript), shQuote(script), shQuote(log_file))
  status <- system2("bash", c("-c", shQuote(cmd)))
  minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  timings <- rbind(timings, data.frame(step = step, minutes = minutes,
                                       status = ifelse(status == 0, "ok", "FAILED")))
  message("==== ", step, " finished in ", minutes, " min (exit status ", status, ")")

  if(status != 0){
    print(timings, row.names = FALSE)
    stop("Step ", step, " failed - see ", log_file, call. = FALSE)
  }

}

message("\nAll requested steps finished.")
print(timings, row.names = FALSE)
