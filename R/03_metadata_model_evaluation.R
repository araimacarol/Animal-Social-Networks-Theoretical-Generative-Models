# =============================================================================
# 03_metadata_model_evaluation.R
# -----------------------------------------------------------------------------
# Evaluate both metadata models (with / without species) against a majority-class
# baseline, reporting macro-F1, balanced accuracy and class-wise precision/recall
# (Reason: accuracy alone might be uninformative when one class dominates).
#
# Uses the RECOVERED, version-portable boosters. If you need to re-export them,
# run recover_booster() ONCE in the recovery environment (older xgboost) -- see
# README > "Model recovery".
#
# Inputs : models/booster_portable*.json, models/preproc_trained_meta*.rds,
#          models/metadata_test.rds, models/metadata_test2.rds
# Outputs: tables/table_meta1.tex, tables/table_meta2.tex
#          outputs/metrics_meta1.csv, outputs/metrics_meta2.csv
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

## ---- One-time recovery (run in the recovery env, then comment out) ----------
# recover_booster(here::here("models", "xgb.final.metadata.fit.new.rds"),
#                 here::here("models", "booster_portable.json"),
#                 here::here("models", "preproc_trained_meta1.rds"))
# recover_booster(here::here("models", "xgb.final.metadata.fit.new2.rds"),
#                 here::here("models", "booster_portable_meta2.json"),
#                 here::here("models", "preproc_trained_meta2.rds"))

evaluate_metadata_model <- function(test_rds, preproc_rds, booster_json, named_json,
                                    caption, label1, label2, tex_out, csv_out) {
  new_data <- readRDS(test_rds)
  preproc  <- readRDS(preproc_rds)

  setup <- finalize_booster(booster_json, preproc, new_data[1, ], named_json)
  res   <- predict_with_labels(new_data, setup$booster, preproc,
                               setup$correct_order, CLASS_ORDER)
  ev    <- evaluate_model(res)

  print(ev$conf_mat)
  print(ev$metrics_table)

  utils::write.csv(ev$metrics_table, csv_out, row.names = FALSE)
  make_latex_table(ev$metrics_table, caption, label1, label2, file = tex_out)
  ev
}

## ---- Model 1: with species --------------------------------------------------
eval_1 <- evaluate_metadata_model(
  test_rds     = here::here("models", "metadata_test.rds"),
  preproc_rds  = here::here("models", "preproc_trained_meta1.rds"),
  booster_json = here::here("models", "booster_portable.json"),
  named_json   = here::here("models", "booster_portable_named.json"),
  caption      = paste("Comprehensive Evaluation Metrics for the Recovered Metadata",
                       "XGBoost Model Against Majority-Class Baseline."),
  label1 = "tab:model_performance", label2 = "tab:species_model_extended_metrics",
  tex_out = here::here("tables",  "table_meta1.tex"),
  csv_out = here::here("outputs", "metrics_meta1.csv"))

## ---- Model 2: without species -----------------------------------------------
eval_2 <- evaluate_metadata_model(
  test_rds     = here::here("models", "metadata_test2.rds"),
  preproc_rds  = here::here("models", "preproc_trained_meta2.rds"),
  booster_json = here::here("models", "booster_portable_meta2.json"),
  named_json   = here::here("models", "booster_portable_named_meta2.json"),
  caption      = paste("Comprehensive Evaluation Metrics for the Recovered Metadata",
                       "XGBoost Model (No Species Feature) Against Majority-Class Baseline."),
  label1 = "tab:model_performance_meta2", label2 = "tab:no_species_model_extended_metrics",
  tex_out = here::here("tables",  "table_meta2.tex"),
  csv_out = here::here("outputs", "metrics_meta2.csv"))

## Class-wise detail behind the supplementary confusion matrices
utils::write.csv(eval_1$model_per_class,
                 here::here("outputs", "classwise_meta1.csv"), row.names = FALSE)
utils::write.csv(eval_2$model_per_class,
                 here::here("outputs", "classwise_meta2.csv"), row.names = FALSE)

message("03_metadata_model_evaluation.R complete.")
