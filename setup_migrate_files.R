# =============================================================================
# setup_migrate_files.R  --  ONE-TIME: copy your existing files into the repo
# -----------------------------------------------------------------------------
# This repo expects files sorted into
# data/ and models/. Run this ONCE to copy them across.
#
# HOW TO USE
#   1. Open Animal-Social-Networks-Theoretical-Generative-Models.Rproj in RStudio.
#   2. Edit OLD_DIR below to point at your old flat working directory.
#   3. source("setup_migrate_files.R")
#   4. Read the report it prints. Anything listed as MISSING must be found
#      before the pipeline will run.
#
# =============================================================================

## ---------------------------------------------------------
OLD_DIR <- "C:/Users/araim/Desktop/ds/R/Workflow/igraphEpi-New"
## -----------------------------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!dir.exists(OLD_DIR)) stop("OLD_DIR does not exist: ", OLD_DIR)

for (d in c("data", "data/edges", "data/Networks-Threshold-Test",
            "models", "outputs", "tables", "figures")) {
  dir.create(here::here(d), recursive = TRUE, showWarnings = FALSE)
}

## Which file goes where. LHS = filename in your old folder, RHS = destination.
manifest <- list(
  ## ---- data/ : raw analysis inputs ----
  data = c(
    "GraphFeatOnAllAnimNets.csv",              # graph features, all empirical nets
    "Network_repository_metaData_combined.csv",# network metadata
    "Badger_Groups.csv"                        # prescribed badger groups
  ),
  ## ---- models/ : model objects + splits ----
  models = c(
    ## RECOVERED, version-portable boosters (these fix the XGBoost reload error)
    "booster_portable_mainmodel.json",
    "preproc_trained_mainmodel.rds",
    "booster_portable.json",                   # metadata model 1 (with species)
    "preproc_trained_meta1.rds",
    "booster_portable_meta2.json",             # metadata model 2 (no species)
    "preproc_trained_meta2.rds",
    ## original fits (still needed by SHAP / hstats in step 07)
    "xgb_final_trained_model.rds",
    "xgb.final.metadata.fit.new.rds",
    "xgb.final.metadata.fit.new2.rds",
    ## prepped recipes + splits
    "metadata.prep.new.rds",
    "metadata.prep.new2.rds",
    "metadata_train.rds",  "metadata_test.rds",
    "metadata_train2.rds", "metadata_test2.rds",
    "df_final_meta_data.rds"
  )
)

copied <- character(); missing <- character()

for (dest in names(manifest)) {
  for (f in manifest[[dest]]) {
    src <- file.path(OLD_DIR, f)
    if (file.exists(src)) {
      file.copy(src, here::here(dest, f), overwrite = TRUE)
      copied <- c(copied, paste0(dest, "/", f))
    } else {
      missing <- c(missing, f)
    }
  }
}

## ---- .edges networks --------------------------------------------------------
## Everything except the threshold-test set goes to data/Empirical-Network-Edges/.
## Original file names are preserved (some carry a ".csv .edges" suffix); the
## pipeline resolves those via resolve_edge_file() rather than renaming them.
edge_dest <- here::here("data", "Empirical-Network-Edges")
if (!dir.exists(edge_dest)) dir.create(edge_dest, recursive = TRUE)
edge_src <- list.files(OLD_DIR, pattern = "\\.edges$", full.names = TRUE)
if (length(edge_src) > 0) {
  file.copy(edge_src, edge_dest, overwrite = TRUE)
  copied <- c(copied, paste0("data/Empirical-Network-Edges/ (", length(edge_src), " .edges files)"))
}

## The threshold-test networks lived in their own subfolder in the old project.
thr_src <- file.path(OLD_DIR, "Networks-Threshold-Test")
if (dir.exists(thr_src)) {
  f <- list.files(thr_src, pattern = "\\.edges$", full.names = TRUE)
  file.copy(f, here::here("data", "Networks-Threshold-Test"), overwrite = TRUE)
  copied <- c(copied, paste0("data/Networks-Threshold-Test/ (", length(f), " files)"))
} else {
  missing <- c(missing, "Networks-Threshold-Test/ (folder of .edges for step 06)")
}

## ---- Report -----------------------------------------------------------------
cat("\n=========================================================\n")
cat("COPIED (", length(copied), "):\n", sep = "")
cat(paste0("  OK  ", copied, collapse = "\n"), "\n")

if (length(missing) > 0) {
  cat("\nMISSING (", length(missing), ") -- find these before running the pipeline:\n", sep = "")
  cat(paste0("  !!  ", missing, collapse = "\n"), "\n")
  cat("\nIf a name differs in your old folder, either rename it to match, or edit\n")
  cat("the `manifest` above. If a RECOVERED booster (*.json / preproc_trained_*)\n")
  cat("is missing, you still need to export it -- see README > 'Model recovery'.\n")
} else {
  cat("\nAll expected files present.\n")
}
cat("=========================================================\n")
cat("\nNext:  source(\"R/02_first_stage_and_metadata.R\")\n\n")
