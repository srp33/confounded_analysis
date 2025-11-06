# Real Data Analysis Script
# Enhanced with functional programming and improved error handling
# Preserves all scientific logic while improving maintainability

options(warn = -1)
suppressMessages(suppressWarnings({
  rm(list=ls())
}))

# Load data and dependencies
load("~/confounded_analysis/scripts/evaluations/robustifying/data/TB_real_data.RData")
source("~/confounded_analysis/scripts/evaluations/robustifying/code/helper.R")
source("~/confounded_analysis/scripts/adjust/gmm_adjust.R")
source("scripts/evaluations/book_chapter/scripts/common_functions.R")
set.seed(123)

# Parse command line arguments
command_args <- commandArgs(trailingOnly=TRUE)

# Determine which analyses to run
if(length(command_args) == 0) {
  cat("Usage: Rscript real_data_analysis.R [3|4|5|6|all] [debug]\n")
  cat("Examples:\n")
  cat("  Rscript real_data_analysis.R 3        # Run 3-study analysis\n")
  cat("  Rscript real_data_analysis.R all      # Run all analyses (3,4,5,6 studies)\n")
  cat("  Rscript real_data_analysis.R 4 debug  # Run 4-study analysis in debug mode\n")
  quit(status=1)
}

analysis_type <- command_args[1]
debug_mode <- length(command_args) > 1 && command_args[2] == "debug"

# Determine which study counts to run
study_counts <- case_when(
  analysis_type == "all" ~ list(c(3, 4, 5, 6)),
  analysis_type %in% c("3", "4", "5", "6") ~ list(as.numeric(analysis_type)),
  TRUE ~ {
    cat("Error: Invalid analysis type. Use 3, 4, 5, 6, or 'all'\n")
    quit(status=1)
  }
)[[1]]

if(analysis_type == "all") {
  cat("Running all analyses: 3, 4, 5, and 6 studies\n")
} else {
  cat(sprintf("Running %d-study analysis\n", study_counts))
}

# Run the specified analyses using purrr::map
analysis_results <- map(study_counts, ~{
  n_studies <- .x
  cat(sprintf("\n=== Starting %d-study analysis ===\n", n_studies))
  start_time <- Sys.time()
  
  result <- tryCatch({
    study_results <- run_analysis_pipeline(
      n_studies = n_studies, 
      dat_lst = dat_lst, 
      label_lst = label_lst, 
      debug_mode = debug_mode
    )
    
    end_time <- Sys.time()
    duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
    
    cat(sprintf("=== Completed %d-study analysis in %.2f minutes ===\n", n_studies, duration))
    
    list(
      n_studies = n_studies,
      success = TRUE,
      duration_mins = duration,
      study_results = study_results,
      error = NULL
    )
    
  }, error = function(e) {
    cat(sprintf("Error in %d-study analysis: %s\n", n_studies, e$message))
    list(
      n_studies = n_studies,
      success = FALSE,
      duration_mins = NA,
      study_results = NULL,
      error = e$message
    )
  })
  
  result
})

names(analysis_results) <- paste0("studies_", study_counts)

# Summary report using functional programming
cat("\n=== ANALYSIS SUMMARY ===\n")

# Count successful analyses
successful_analyses <- sum(map_lgl(analysis_results, "success"))
total_analyses <- length(analysis_results)

cat(sprintf("Successful analyses: %d/%d\n", successful_analyses, total_analyses))

# Report timing for successful analyses
successful_results <- keep(analysis_results, "success")
if(length(successful_results) > 0) {
  total_time <- sum(map_dbl(successful_results, "duration_mins"), na.rm = TRUE)
  cat(sprintf("Total runtime: %.2f minutes\n", total_time))
  
  # Individual timings
  walk(successful_results, ~{
    cat(sprintf("  %d studies: %.2f minutes\n", .x$n_studies, .x$duration_mins))
  })
}

# Report errors if any
failed_results <- keep(analysis_results, ~!.x$success)
if(length(failed_results) > 0) {
  cat("\nErrors encountered:\n")
  walk(failed_results, ~{
    cat(sprintf("  %d studies: %s\n", .x$n_studies, .x$error))
  })
}

cat("\nAll requested analyses completed.\n")

# Optional: Save results summary
if(!debug_mode && successful_analyses > 0) {
  results_summary <- map_dfr(analysis_results, ~{
    tibble(
      n_studies = .x$n_studies,
      success = .x$success,
      duration_mins = .x$duration_mins %||% NA,
      error = .x$error %||% NA
    )
  })
  
  summary_file <- sprintf("~/confounded_analysis/scripts/evaluations/robustifying/analysis_summary_%s.csv", 
                         format(Sys.time(), "%Y%m%d_%H%M%S"))
  write.csv(results_summary, summary_file, row.names = FALSE)
  cat(sprintf("Results summary saved to: %s\n", summary_file))
}