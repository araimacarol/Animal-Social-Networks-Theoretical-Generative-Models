# =============================================================================
# run_all.R  --  master script: runs the whole pipeline end to end
# -----------------------------------------------------------------------------
# Usage:
#
#   Rscript run_all.R          # run everything
#   Rscript run_all.R 05 06    # run only steps 05 and 06
#
# Or from RStudio: open the .Rproj, then source("run_all.R").
#
# Working directory does not matter: {here} anchors every path to the project
# root (the folder containing the .Rproj / .git file), so this works whether you
# launch from the root, from R/, or from anywhere inside the project.
#
# Steps 04 (grouped CV) and 06 (threshold sensitivity) are the slow ones -- they
# refit / re-score models many times.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

steps <- c(
  "02" = "02_first_stage_and_metadata.R",
  "03" = "03_metadata_model_evaluation.R",
  "04" = "04_grouped_cross_validation.R",
  "05" = "05_group_level_community_detection.R",
  "06" = "06_threshold_sensitivity.R",
  "07" = "07_interpretation_shap_hstats.R"
)

args <- commandArgs(trailingOnly = TRUE)
to_run <- if (length(args) == 0) names(steps) else intersect(args, names(steps))
if (length(to_run) == 0) stop("No valid steps. Choose from: ", paste(names(steps), collapse = ", "))

for (s in to_run) {
  message("\n==============================================================")
  message(">>> STEP ", s, ": ", steps[[s]])
  message("==============================================================")
  source(here::here("R", steps[[s]]), echo = FALSE)
}

message("\nPipeline finished. All generated tables were written to outputs/.")
