# =============================================================================
# 04_grouped_cross_validation.R
# -----------------------------------------------------------------------------
# Leave-one-group-out cross-validation for the metadata models.
#
# A random 70:30 split overestimates performance because many
# networks are not independent (same species, study, season, colony, protocol).
# We hold out an entire group at a time -- each taxonomic class, interaction
# type, data-collection method and captivity status exhaustively, plus the
# majority species for the species-included model -- and re-evaluate.
#
# Inputs : models/df_final_meta_data.rds
# Outputs: outputs/pipeline_overall_performance2.csv
#          outputs/pipeline_classwise_performance2.csv  (+ .rds copies)
#
# NOTE: this refits and re-tunes XGBoost once per held-out group. Expect it to
# run for a while; it is the slowest step in the pipeline.
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

pipeline_results <- run_ecological_network_pipeline(
  data_path        = here::here("models", "df_final_meta_data.rds"),
  seed             = 123,
  top_n_species    = 5,    # species has too many levels for exhaustive LOGO
  tuning_grid_size = 15
)

cat("\n### OVERALL PERFORMANCE (random vs grouped validation) ###\n")
print(knitr::kable(pipeline_results$overall_performance, format = "markdown", digits = 4))

cat("\n### CLASS-WISE METRICS ###\n")
print(knitr::kable(pipeline_results$classwise_performance, format = "markdown", digits = 4))

saveRDS(pipeline_results$overall_performance,
        here::here("outputs", "pipeline_overall_performance2.rds"))
saveRDS(pipeline_results$classwise_performance,
        here::here("outputs", "pipeline_classwise_performance2.rds"))

utils::write.csv(pipeline_results$overall_performance,
                 here::here("outputs", "pipeline_overall_performance2.csv"), row.names = FALSE)
utils::write.csv(pipeline_results$classwise_performance,
                 here::here("outputs", "pipeline_classwise_performance2.csv"), row.names = FALSE)

message("04_grouped_cross_validation.R complete.")
