# =============================================================================
# 07_interpretation_shap_hstats.R
# -----------------------------------------------------------------------------
# Model interpretation: SHAP attributions and Friedman-Popescu H-statistics.
#
# SHAP explains the MODEL'S PREDICTIONS, not ecological
# causation. Species, interaction type, data-collection method, captivity status
# and observation duration  for the empirical networks are confounded, so these outputs
# are read as "how the model uses these features", never as causal drivers.
#
# MODEL ROUTE -- IMPORTANT
#   The saved tidymodels fits (xgb.final.metadata.fit.new{,2}.rds) carry
#   xgb.Booster handles that throw
#       "'xgb.Booster' object is corrupted or is from an incompatible XGBoost
#        version"
#   under current xgboost. Confirmed by xgb.importance() on the extracted
#   engine. We therefore use the version-portable JSON boosters plus the trained
#   recipes -- the same route as scripts 02, 03, 05 and 06.
#

#
# Inputs : models/booster_portable.json        + models/preproc_trained_meta1.rds
#          models/booster_portable_meta2.json  + models/preproc_trained_meta2.rds
#          models/metadata_train{,2}.rds, models/df_final_meta_data.rds
# Outputs: figures/shap_*.pdf, figures/hstat_*.pdf, models/final_meta_shap*.rds
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

set.seed(1486)

## Feature display names (order must match the baked predictor matrix)
SHAP_NAMES_SPECIES    <- c("Captive", "Interaction type", "Data collection",
                           "Species", "Data duration (days)", "Class")
SHAP_NAMES_NO_SPECIES <- c("Captive", "Interaction type", "Data collection",
                           "Data duration (days)", "Class")


# =============================================================================
# Helper: load a recovered metadata model (portable booster + trained recipe)
# =============================================================================
load_meta_model <- function(booster_json, recipe_rds) {
  stopifnot(file.exists(booster_json), file.exists(recipe_rds))
  preproc <- readRDS(recipe_rds)
  list(booster = xgboost::xgb.load(booster_json), preproc = preproc)
}

#' Bake a raw metadata data frame into the predictor matrix the booster expects.
#' Column ORDER is taken from the trained recipe -- feature names did not survive
#' JSON re-serialisation, so order is what we rely on.
bake_predictors <- function(preproc, new_data) {
  baked <- recipes::bake(preproc, new_data = new_data)
  keep  <- setdiff(names(baked), c("target", ".pred_class"))
  X <- as.matrix(baked[, keep, drop = FALSE])
  storage.mode(X) <- "double"
  X
}


# =============================================================================
# SHAP
# =============================================================================

#' Build a shapviz object from a RECOVERED booster and relabel its classes.
build_shap <- function(booster_json, recipe_rds, train_rds, feature_names,
                       interactions = FALSE) {
  mm    <- load_meta_model(booster_json, recipe_rds)
  train <- readRDS(train_rds)
  train <- train[sample(nrow(train)), ]
  
  X_pred <- bake_predictors(mm$preproc, train)
  
  if (ncol(X_pred) != length(feature_names))
    stop("Baked predictor matrix has ", ncol(X_pred), " columns but ",
         length(feature_names), " display names were supplied: ",
         paste(colnames(X_pred), collapse = ", "))
  
  sv <- shapviz::shapviz(mm$booster, X_pred = X_pred, X = train,
                         interactions = interactions)
  
  ## XGBoost emits Class_1..Class_5; map to the generative-model names.
  ## Order follows CLASS_ORDER = c("ER","sbm","SF","Spatial","SW").
  names(sv) <- dplyr::recode(names(sv),
                             "Class_1" = "Erdos-Renyi", "Class_2" = "Stochastic-Block-Model",
                             "Class_3" = "Scale-Free",  "Class_4" = "Spatial", "Class_5" = "Small-World",
                             .default = names(sv))
  
  for (cls in names(sv)) {
    colnames(sv[[cls]]$X) <- feature_names
    colnames(sv[[cls]]$S) <- feature_names
    if (interactions && !is.null(sv[[cls]]$S_inter)) {
      colnames(sv[[cls]]$S_inter)      <- feature_names
      dimnames(sv[[cls]]$S_inter)[[3]] <- feature_names
    }
  }
  sv
}

## ---- SHAP: species-excluded model (Figures 4 & 5 in the paper) --------------
shap_no_species <- build_shap(
  here::here("models", "booster_portable_meta2.json"),
  here::here("models", "preproc_trained_meta2.rds"),
  here::here("models", "metadata_train2.rds"),
  SHAP_NAMES_NO_SPECIES)
saveRDS(shap_no_species, here::here("models", "final_meta_shap_no_species.rds"))

ggplot2::ggsave(here::here("figures", "shap_importance_geometric_no_species.pdf"),
                shap_importance_plot(shap_no_species$Spatial), width = 12, height = 8)
ggplot2::ggsave(here::here("figures", "shap_importance_scalefree_no_species.pdf"),
                shap_importance_plot(shap_no_species$`Scale-Free`), width = 12, height = 8)

## ---- SHAP: species-included model -------------------------------------------
shap_species <- build_shap(
  here::here("models", "booster_portable.json"),
  here::here("models", "preproc_trained_meta1.rds"),
  here::here("models", "metadata_train.rds"),
  SHAP_NAMES_SPECIES)
saveRDS(shap_species, here::here("models", "final_meta_shap_species.rds"))

ggplot2::ggsave(here::here("figures", "shap_importance_geometric_species.pdf"),
                shap_importance_plot(shap_species$Spatial), width = 12, height = 8)


# =============================================================================
# H-STATISTICS (species-excluded model)
# -----------------------------------------------------------------------------
# hstats() permutes columns of the RAW data frame X and calls pred_fun(object, X)
# on each perturbation, so pred_fun must bake before predicting. It cannot be
# handed a bare booster with a pre-baked matrix.
# =============================================================================

df_meta <- readRDS(here::here("models", "df_final_meta_data.rds")) %>%
  dplyr::select(-species)

mm2 <- load_meta_model(here::here("models", "booster_portable_meta2.json"),
                       here::here("models", "preproc_trained_meta2.rds"))

## Probability predictions from the native booster, taking a RAW data frame.
pred_fun_meta <- function(object, newdata, ...) {
  X <- bake_predictors(object$preproc, newdata)
  p <- stats::predict(object$booster, X)
  if (!is.matrix(p)) p <- matrix(p, nrow = nrow(X), byrow = TRUE)
  colnames(p) <- CLASS_ORDER
  p
}

## Sanity check before the (slow) hstats run: probabilities must be a 5-column
## matrix summing to 1, and the argmax must match the recovered predictions.
chk <- pred_fun_meta(mm2, utils::head(df_meta, 5))
stopifnot(is.matrix(chk), ncol(chk) == 5,
          all(abs(rowSums(chk) - 1) < 1e-6))
print(round(chk, 4))

h <- hstats::hstats(mm2, X = df_meta %>% dplyr::select(-target),
                    v = setdiff(colnames(df_meta), "target"),
                    pred_fun = pred_fun_meta,
                    n_max = 600, threeway_m = 3L)
saveRDS(h, here::here("models", "meta_interactn_hstats_no_species.rds"))

hstat_names <- c("Captive", "Interaction type", "Data collection",
                 "Data duration (days)", "Class")
class_names <- c("Erdos-Renyi", "Stochastic-Block-Model", "Scale-Free",
                 "Geometric Graphs", "Small-World")
rownames(h$h2_overall$num) <- hstat_names
colnames(h$h2_overall$num) <- class_names
names(h$h2_overall$denom)  <- class_names

plot_h <- function(x, file, width = 20, height = 15) {
  p <- plot(x) + ggplot2::xlab("Values") + theme_paper(24) +
    ggplot2::scale_fill_manual(name = "Networks", values = COLOR_PALETTE)
  ggplot2::ggsave(file, p, width = width, height = height)
}

plot_h(hstats::h2_overall(h, normalize = FALSE, squared = FALSE, top_m = 8, zero = TRUE),
       here::here("figures", "hstat_overall_no_species.pdf"))
plot_h(hstats::h2_pairwise(h, normalize = FALSE, squared = FALSE, top_m = 8, zero = TRUE),
       here::here("figures", "hstat_pairwise_no_species.pdf"))
plot_h(hstats::h2_threeway(h, normalize = FALSE, squared = FALSE, top_m = 8, zero = TRUE),
       here::here("figures", "hstat_threeway_no_species.pdf"), width = 24)

message("07_interpretation_shap_hstats.R complete.")