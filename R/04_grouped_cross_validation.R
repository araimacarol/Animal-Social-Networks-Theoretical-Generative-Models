# =============================================================================
# 04_grouped_cross_validation.R
# -----------------------------------------------------------------------------
# Random-split and grouped (leave-one-group-out) validation of the two
# metadata models: Model 1 (with species) and Model 2 (without species).
#
#   1. Identical splits. Both models are fitted and tested on exactly the same
#      rows in every comparison: one stratified 70:30 random split, and the same
#      held-out groups.
#   2. Leave-one-species-out for BOTH models (five most frequent species),
#      reported separately from hold-outs of an entire metadata category
#      (taxonomic class, interaction type, data-collection method, captivity).
#   3. Majority-class baseline taken from each TRAINING fold (the class most
#      common in training, scored on the test fold).
#   4. Unseen levels handled explicitly. Factor levels absent from a training
#      fold are dropped before fitting (droplevels), and step_novel() maps any
#      such level in the test fold to a single "new" level. Previously a
#      held-out level kept its integer code from the full factor, so test rows
#      carried a code no training row had. The share of test rows with a level
#      unseen in training is reported for each predictor.
#   The model specification is otherwise unchanged (XGBoost, step_integer,
#   tune_race_anova over a 15-point grid, accuracy as the selection metric).
#
# Inputs : models/df_final_meta_data.rds
# Outputs: outputs/point4_cv_overall.csv
#          outputs/point4_cv_classwise.csv
#
# NOTE: refits and re-tunes XGBoost twice per fold. Slow.
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

SEED          <- 123
TOP_N_SPECIES <- 5
GRID_SIZE     <- 15
CATEGORICAL   <- c("species", "class", "interaction_type", "data_collection", "captive")

df <- readRDS(here::here("models", "df_final_meta_data.rds")) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(c("target", CATEGORICAL)), as.factor))

MODELS <- list(
  With_Species    = setdiff(names(df), "target"),
  Without_Species = setdiff(names(df), c("target", "species")))

xgb_spec <- parsnip::boost_tree(
  mtry = tune::tune(), trees = tune::tune(),
  learn_rate = tune::tune(), tree_depth = tune::tune()) %>%
  parsnip::set_engine("xgboost") %>%
  parsnip::set_mode("classification")

## ---- Metrics (same per-class convention as before: undefined -> 0) ----------
class_metrics <- function(truth, pred, lev) {
  purrr::map_df(lev, function(cl) {
    tp <- sum(pred == cl & truth == cl); np <- sum(pred == cl); na <- sum(truth == cl)
    p <- if (np == 0) 0 else tp / np; r <- if (na == 0) 0 else tp / na
    tibble::tibble(Target_Class = cl, Class_Precision = p, Class_Recall = r,
                   Class_F1 = if (p + r == 0) 0 else 2 * p * r / (p + r), Class_Support = na)
  })
}

## ---- Fit one model on one fold and score it ---------------------------------
fit_and_score <- function(train, test, predictors, model_id, validation, held_out, fold_seed) {
  cat(sprintf("   --> %s | %s | held out: %s\n", model_id, validation, held_out))
  train <- droplevels(train[, c("target", predictors)])
  test  <- test[, c("target", predictors)]

  rec <- recipes::recipe(target ~ ., data = train) %>%
    recipes::step_novel(recipes::all_nominal_predictors()) %>%
    recipes::step_integer(recipes::all_nominal_predictors())
  wf  <- workflows::workflow() %>% workflows::add_model(xgb_spec) %>% workflows::add_recipe(rec)

  set.seed(fold_seed)                         # same tuning seed for both models
  v     <- max(2, min(5, floor(nrow(train) / 10)))
  folds <- rsample::vfold_cv(train, v = v, strata = target)
  tuned <- tryCatch(
    finetune::tune_race_anova(wf, resamples = folds, grid = GRID_SIZE,
      control = finetune::control_race(verbose_elim = FALSE),
      metrics = yardstick::metric_set(yardstick::accuracy)),
    error = function(e) { message("   tuning failed: ", conditionMessage(e)); NULL })
  if (is.null(tuned)) return(NULL)
  fit  <- parsnip::fit(tune::finalize_workflow(wf, tune::select_best(tuned, metric = "accuracy")), train)
  pred <- as.character(stats::predict(fit, test)$.pred_class)

  lev   <- levels(df$target)
  truth <- as.character(test$target)
  maj   <- names(which.max(table(train$target)))          # training-fold majority
  cw    <- class_metrics(truth, pred, lev)
  cw_b  <- class_metrics(truth, rep(maj, length(truth)), lev)

  unseen <- vapply(intersect(predictors, CATEGORICAL), function(v)
    round(100 * mean(!as.character(test[[v]]) %in% levels(train[[v]])), 1), numeric(1))

  list(
    overall = tibble::tibble(
      Model = model_id, Validation = validation, Held_Out_Group = held_out,
      N_Train = nrow(train), N_Test = nrow(test),
      Accuracy = round(mean(pred == truth), 4),
      Macro_F1 = round(mean(cw$Class_F1), 4),
      Baseline_Class = maj,
      Baseline_Accuracy = round(mean(truth == maj), 4),
      Baseline_Macro_F1 = round(mean(cw_b$Class_F1), 4),
      N_Distinct_Predicted = dplyr::n_distinct(pred),
      Pct_Test_Unseen = paste(sprintf("%s=%.1f", names(unseen), unseen), collapse = "; ")),
    classwise = dplyr::mutate(cw, Model = model_id, Validation = validation,
                              Held_Out_Group = held_out, .before = 1))
}

## ---- Define every fold once, then fit both models on it ---------------------
set.seed(SEED)
split <- rsample::initial_split(df, prop = 0.7, strata = target)
folds_def <- list(list(validation = "Random_70_30_Split", held_out = "None",
                       train = rsample::training(split), test = rsample::testing(split)))

top_species <- names(sort(table(df$species), decreasing = TRUE))[seq_len(TOP_N_SPECIES)]
for (sp in top_species)
  folds_def[[length(folds_def) + 1]] <- list(validation = "Leave_One_Species_Out", held_out = sp,
    train = df[df$species != sp, ], test = df[df$species == sp, ])

for (col in setdiff(CATEGORICAL, "species")) {
  for (grp in names(sort(table(df[[col]]), decreasing = TRUE))) {
    tr <- df[df[[col]] != grp, ]; te <- df[df[[col]] == grp, ]
    if (nrow(te) == 0 || dplyr::n_distinct(tr$target) < 2) next
    folds_def[[length(folds_def) + 1]] <- list(
      validation = paste0("Leave_One_", col, "_Out"), held_out = grp, train = tr, test = te)
  }
}

overall <- list(); classwise <- list()
for (i in seq_along(folds_def)) {
  f <- folds_def[[i]]
  for (m in names(MODELS)) {
    res <- fit_and_score(f$train, f$test, MODELS[[m]], m, f$validation, f$held_out,
                         fold_seed = SEED + i)
    if (!is.null(res)) { overall[[length(overall) + 1]] <- res$overall
                         classwise[[length(classwise) + 1]] <- res$classwise }
  }
}
overall   <- dplyr::bind_rows(overall) %>%
  dplyr::mutate(Scheme = dplyr::case_when(
    Validation == "Random_70_30_Split"    ~ "Random split",
    Validation == "Leave_One_Species_Out" ~ "Species hold-out",
    TRUE                                  ~ "Metadata-category hold-out"), .before = 1)
classwise <- dplyr::bind_rows(classwise)

print(as.data.frame(overall), row.names = FALSE)

## ---- Side-by-side comparison of the two models on identical folds -----------
paired <- overall %>%
  dplyr::select(Scheme, Validation, Held_Out_Group, N_Test, Baseline_Class, Baseline_Accuracy,
                Model, Accuracy, Macro_F1) %>%
  tidyr::pivot_wider(names_from = Model, values_from = c(Accuracy, Macro_F1)) %>%
  dplyr::mutate(Diff_Accuracy = Accuracy_With_Species - Accuracy_Without_Species)
cat("\nPaired comparison (identical folds):\n"); print(as.data.frame(paired), row.names = FALSE)

utils::write.csv(overall,   here::here("outputs", "point4_cv_overall.csv"),   row.names = FALSE)
utils::write.csv(classwise, here::here("outputs", "point4_cv_classwise.csv"), row.names = FALSE)
utils::write.csv(paired,    here::here("outputs", "point4_cv_paired.csv"),    row.names = FALSE)

message("04_grouped_cross_validation.R complete.")
