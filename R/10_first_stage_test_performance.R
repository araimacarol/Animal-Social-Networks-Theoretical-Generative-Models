# =============================================================================
# 10_first_stage_test_performance.R
# -----------------------------------------------------------------------------
# Main Table 2: performance of the first-stage classifier on the held-out
# simulated test set (132,643 networks), computed with the RECOVERED booster
# (models/booster_portable_mainmodel.json + preproc_trained_mainmodel.rds).
#
# Input : TEST_RDS -- the held-out simulated test networks (feature columns as in
#         training plus the true class). Set the path below to that file; it is
#         the 30% test partition of the 442,142 simulated networks.
# Output: outputs/table2_first_stage_confusion.csv  (rows = predicted, cols = true)
#         outputs/table2_first_stage_metrics.csv    (class-wise precision/recall,
#                                                   accuracy, Hand-Till AUC, log-loss)
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

TEST_RDS <- here::here("models", "df_test_new.rds")   # <- set to the test-set file
if (!file.exists(TEST_RDS)) stop("Test-set file not found: ", TEST_RDS)

model <- load_main_model()
test  <- readRDS(TEST_RDS)
lab   <- intersect(c("graph_name", "target", "GraphName", "class"), names(test))[1]
truth <- factor(as.character(test[[lab]]), levels = CLASS_ORDER)
stopifnot(!anyNA(truth))
cat("Test networks:", nrow(test), "(expect 132,643)\n")

## Class probabilities from the recovered booster (same route as predict_main_model)
feats <- test %>% dplyr::mutate(dplyr::across(dplyr::where(is.integer), as.double))
baked <- recipes::bake(model$preproc, new_data = feats)
keep  <- setdiff(names(baked), c("graph_name", "target", "class", ".pred_class", "GraphName"))
X <- as.matrix(baked[, keep, drop = FALSE]); storage.mode(X) <- "double"
P <- stats::predict(model$booster, X)
if (!is.matrix(P)) P <- matrix(P, nrow = nrow(X), byrow = TRUE)
colnames(P) <- CLASS_ORDER
pred <- factor(CLASS_ORDER[max.col(P, ties.method = "first")], levels = CLASS_ORDER)

cm <- table(Predicted = pred, Truth = truth)
print(cm)

res <- tibble::tibble(truth = truth, predicted_class = pred) %>% dplyr::bind_cols(tibble::as_tibble(P))
metrics <- tibble::tibble(
  class     = CLASS_ORDER,
  precision = round(diag(cm) / rowSums(cm), 4),
  recall    = round(diag(cm) / colSums(cm), 4),
  support   = as.integer(colSums(cm)))
overall <- tibble::tibble(
  accuracy      = round(mean(pred == truth), 4),
  roc_auc_hand_till = round(yardstick::roc_auc(res, truth = truth, dplyr::all_of(CLASS_ORDER),
                                                estimator = "hand_till")$.estimate, 4),
  log_loss      = signif(yardstick::mn_log_loss(res, truth = truth, dplyr::all_of(CLASS_ORDER))$.estimate, 3))
print(metrics); print(overall)

utils::write.csv(as.data.frame.matrix(cm), here::here("outputs", "table2_first_stage_confusion.csv"))
utils::write.csv(dplyr::bind_cols(metrics, overall[rep(1, nrow(metrics)), ]),
                 here::here("outputs", "table2_first_stage_metrics.csv"), row.names = FALSE)

message("10_first_stage_test_performance.R complete.")
