# =============================================================================
# run_all.R  --  master script: runs the whole pipeline end to end
# -----------------------------------------------------------------------------
# Usage:
#   Rscript run_all.R          # run everything
#   Rscript run_all.R 05 06    # run only steps 05 and 06
# Or from RStudio: open the .Rproj, then source("run_all.R").
#
# {here} anchors every path to the project root, so the working directory does
# not matter. Every classification uses the RECOVERED models
# (models/booster_portable*.json + models/preproc_trained_*.rds).
#
# Slow steps: 04 (refits both metadata models on every fold), 06 (set RERUN),
# PCA (reads the training data unless point2_export_derived.R has been run).
# Run separately (not here): sbm-sim-different-communities.R (simulations
# behind S16/S17), point2_export_derived.R (one-off, needs the large training
# files), 10 (needs the simulated test-set file).
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

steps <- c(
  "02"  = "02_first_stage_and_metadata.R",       # Fig 3, Fig S2, Tables S2-S5
  "03"  = "03_metadata_model_evaluation.R",      # Tables 5, 6, S7
  "04"  = "04_grouped_cross_validation.R",       # Tables 7, S8
  "05"  = "05_group_level_community_detection.R",# Table S6, Results ranges
  "06"  = "06_threshold_sensitivity.R",          # Tables 3, S1, S23, S24
  "07"  = "07_interpretation_shap_hstats.R",     # Figs 4, 5, S3-S11
  "08"  = "08_sbm_simulation_tables.R",          # Tables S16, S17
  "09"  = "09_group_level_tables.R",             # Tables 4, S12-S14
  "PCA" = "PCA_Applicability.R",                 # Fig 2, Fig S1, Tables S18-S22, S25
  "11"  = "11_illustrative_figures.R"            # Fig 1, Fig S12
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
message("\nPipeline finished. Tables in outputs/, figures in figures/.")
