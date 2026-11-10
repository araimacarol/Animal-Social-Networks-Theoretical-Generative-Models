# =============================================================================
# 07_interpretation_shap_hstats.R
# -----------------------------------------------------------------------------
# Model interpretation: SHAP attributions and Friedman-Popescu H-statistics for
# the metadata models (Model 1 with species, Model 2 without species), and the
# code for every SHAP / H-statistic figure in the paper and supplement.
#
# SHAP explains the MODEL'S PREDICTIONS, not ecological causation. Species,
# interaction type, data-collection method, captivity status and observation
# duration are confounded across the empirical networks, so these outputs are
# read as "how the model uses these features", never as causal drivers.
#
# MODEL ROUTE
#   The saved tidymodels fits (xgb.final.metadata.fit.new{,2}.rds) carry
#   xgb.Booster handles that fail to load under current xgboost ("object is
#   corrupted or is from an incompatible XGBoost version"). We therefore use the
#   version-portable JSON boosters plus the trained recipes -- the recovered
#   models evaluated in Tables 5, 6 and S7 (corrected Model 2), and the same
#   route as scripts 02, 03, 05 and 06. (Manuscript: Tables 5, 6 and S7.)
#   H-statistics use all 636 networks in df_final_meta_data.rds (the data the
#   statistic permutes); SHAP values use the training sets of each model.
#
# WHAT IS PLOTTED 
#   - SHAP values: one value per network (training set), feature and class, on
#     the model's raw output scale, i.e. log-odds (margin) for each class.
#   - Importance plots: bar = mean absolute SHAP value of the feature.
#   - Dependence plots: one point per network at its SHAP value; the black
#     diamond marks the level's mean SHAP value; n = networks per level. Values are not
#     summed, so a level's position does not depend on how common it is.
#   - Waterfall plot: SHAP values of one network, from the base value (mean
#     log-odds) to that network's log-odds for the class.
#   - Two-way plots: one point per network at its pair of feature levels,
#     coloured by the SHAP interaction value of the pair (log-odds).
#   - H-statistics: computed on the predicted-probability scale.
#   - Plotting data written to outputs/ (per-network SHAP values, per-level
#     summaries, importance values).
#
# Inputs : models/booster_portable.json        + models/preproc_trained_meta1.rds
#          models/booster_portable_meta2.json  + models/preproc_trained_meta2.rds
#          models/metadata_train{,2}.rds, models/df_final_meta_data.rds
# Outputs: figures/*.pdf (file names as used in the manuscript),
#          models/final_meta_shap*.rds, models/meta_interactn_hstats_no_species.rds,
#          outputs/shap_values_*.csv, outputs/shap_level_summary.csv,
#          outputs/shap_importance.csv
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

set.seed(1486)

## Feature display names (order must match the baked predictor matrix)
SHAP_NAMES_SPECIES    <- c("Captive", "Interaction type", "Data collection",
                           "Species", "Data duration (days)", "Class")
SHAP_NAMES_NO_SPECIES <- c("Captive", "Interaction type", "Data collection",
                           "Data duration (days)", "Class")
CLASS_LABEL <- c("Spatial" = "geometric graph", "Scale-Free" = "scale-free")

## Network used for the waterfall plot (row of metadata_train2.rds, before
## shuffling). Set to the network shown in the supplementary waterfall figure.
WATERFALL_ROW <- 1

fig <- function(f) here::here("figures", f)

# =============================================================================
# Helpers: load a recovered metadata model and bake predictors
# =============================================================================
load_meta_model <- function(booster_json, recipe_rds) {
  stopifnot(file.exists(booster_json), file.exists(recipe_rds))
  list(booster = xgboost::xgb.load(booster_json), preproc = readRDS(recipe_rds))
}

#' Bake a raw metadata data frame into the predictor matrix the booster expects.
#' Column ORDER is taken from the trained recipe.
bake_predictors <- function(preproc, new_data) {
  baked <- recipes::bake(preproc, new_data = new_data)
  keep  <- setdiff(names(baked), c("target", ".pred_class"))
  X <- as.matrix(baked[, keep, drop = FALSE])
  storage.mode(X) <- "double"
  X
}

# =============================================================================
# 1. SHAP (with interaction values, needed for the two-way plots)
# =============================================================================
build_shap <- function(booster_json, recipe_rds, train_rds, feature_names, shuffle = TRUE) {
  mm    <- load_meta_model(booster_json, recipe_rds)
  train <- readRDS(train_rds)
  if (shuffle) train <- train[sample(nrow(train)), ]
  X_pred <- bake_predictors(mm$preproc, train)
  if (ncol(X_pred) != length(feature_names))
    stop("Baked predictor matrix has ", ncol(X_pred), " columns but ",
         length(feature_names), " display names were supplied: ",
         paste(colnames(X_pred), collapse = ", "))
  X_raw <- as.data.frame(train[, colnames(X_pred)])
  
  sv <- shapviz::shapviz(mm$booster, X_pred = X_pred, X = X_raw, interactions = TRUE)
  ## XGBoost emits Class_1..Class_5, in CLASS_ORDER = c("ER","sbm","SF","Spatial","SW").
  names(sv) <- dplyr::recode(names(sv),
                             "Class_1" = "Erdos-Renyi", "Class_2" = "Stochastic-Block-Model",
                             "Class_3" = "Scale-Free",  "Class_4" = "Spatial", "Class_5" = "Small-World",
                             .default = names(sv))
  for (cls in names(sv)) {
    colnames(sv[[cls]]$X) <- feature_names
    colnames(sv[[cls]]$S) <- feature_names
    if (!is.null(sv[[cls]]$S_inter)) {
      colnames(sv[[cls]]$S_inter)      <- feature_names
      dimnames(sv[[cls]]$S_inter)[[3]] <- feature_names
    }
    ## readable level labels (e.g. "physical_contact" -> "physical contact")
    sv[[cls]]$X[] <- lapply(sv[[cls]]$X, function(v)
      if (is.character(v) || is.factor(v)) factor(gsub("_", " ", as.character(v))) else v)
  }
  sv
}

shap_no_species <- build_shap(here::here("models", "booster_portable_meta2.json"),
                              here::here("models", "preproc_trained_meta2.rds"),
                              here::here("models", "metadata_train2.rds"),
                              SHAP_NAMES_NO_SPECIES)
shap_species    <- build_shap(here::here("models", "booster_portable.json"),
                              here::here("models", "preproc_trained_meta1.rds"),
                              here::here("models", "metadata_train.rds"),
                              SHAP_NAMES_SPECIES)
saveRDS(shap_no_species, here::here("models", "final_meta_shap_no_species.rds"))
saveRDS(shap_species,    here::here("models", "final_meta_shap_species.rds"))

# =============================================================================
# 1b. H-STATISTICS (species-excluded model), probability scale
# Run directly after the two SHAP builds, as before, so the random subsample of
# 600 networks (and hence the H values) is unchanged; plotting below uses the RNG.
# -----------------------------------------------------------------------------
# hstats() permutes columns of the RAW data frame X and calls pred_fun(object, X)
# on each perturbation, so pred_fun must bake before predicting.
# =============================================================================
df_meta <- readRDS(here::here("models", "df_final_meta_data.rds")) %>%
  dplyr::select(-species)
mm2 <- load_meta_model(here::here("models", "booster_portable_meta2.json"),
                       here::here("models", "preproc_trained_meta2.rds"))

pred_fun_meta <- function(object, newdata, ...) {
  X <- bake_predictors(object$preproc, newdata)
  p <- stats::predict(object$booster, X)
  if (!is.matrix(p)) p <- matrix(p, nrow = nrow(X), byrow = TRUE)
  colnames(p) <- CLASS_ORDER
  p
}
chk <- pred_fun_meta(mm2, utils::head(df_meta, 5))
stopifnot(is.matrix(chk), ncol(chk) == 5, all(abs(rowSums(chk) - 1) < 1e-6))

h <- hstats::hstats(mm2, X = df_meta %>% dplyr::select(-target),
                    v = setdiff(colnames(df_meta), "target"),
                    pred_fun = pred_fun_meta, n_max = 600, threeway_m = 3L)
saveRDS(h, here::here("models", "meta_interactn_hstats_no_species.rds"))

hstat_names <- c("Captive", "Interaction type", "Data collection",
                 "Data duration (days)", "Class")
class_names <- c("Erdos-Renyi", "Stochastic-Block-Model", "Scale-Free",
                 "Geometric Graphs", "Small-World")
rownames(h$h2_overall$num) <- hstat_names
colnames(h$h2_overall$num) <- class_names
names(h$h2_overall$denom)  <- class_names

plot_h <- function(x, file, width = 20, height = 15) {
  p <- plot(x) + ggplot2::xlab("H statistic (probability scale)") + theme_paper(24) +
    ggplot2::scale_fill_manual(name = "Networks", values = COLOR_PALETTE)
  ggplot2::ggsave(file, p, width = width, height = height)
}
## file names as used in the supplement
plot_h(hstats::h2_overall(h,  normalize = FALSE, squared = FALSE, top_m = 8, zero = TRUE),
       fig("hstat_overall_nonspecies.pdf"))
plot_h(hstats::h2_pairwise(h, normalize = FALSE, squared = FALSE, top_m = 8, zero = TRUE),
       fig("hstat_pairwise_nonspecies.pdf"))
plot_h(hstats::h2_threeway(h, normalize = FALSE, squared = FALSE, top_m = 8, zero = TRUE),
       fig("hstat_threeway_nonspecies.pdf"), width = 24)


# =============================================================================
# 2. Plotting data (written to outputs/)
# =============================================================================
shap_long <- function(sv_cls, model, cls) {
  S <- shapviz::get_shap_values(sv_cls); X <- shapviz::get_feature_values(sv_cls)
  dplyr::bind_rows(lapply(colnames(S), function(f) tibble::tibble(
    model = model, class = cls, network = seq_len(nrow(S)), feature = f,
    level = as.character(X[[f]]), shap_logodds = S[, f])))
}
importance <- list(); level_summary <- list()
for (spec in list(list(shap_species, "Model1_with_species"),
                  list(shap_no_species, "Model2_without_species"))) {
  for (cls in c("Spatial", "Scale-Free")) {
    d <- shap_long(spec[[1]][[cls]], spec[[2]], cls)
    utils::write.csv(d, here::here("outputs", sprintf("shap_values_%s_%s.csv", spec[[2]], cls)),
                     row.names = FALSE)
    importance[[length(importance) + 1]] <- d %>%
      dplyr::group_by(model, class, feature) %>%
      dplyr::summarise(mean_abs_shap = mean(abs(shap_logodds)), .groups = "drop") %>%
      dplyr::mutate(share_pct = round(100 * mean_abs_shap / sum(mean_abs_shap), 1))
    level_summary[[length(level_summary) + 1]] <- d %>%
      dplyr::filter(feature != "Data duration (days)") %>%
      dplyr::group_by(model, class, feature, level) %>%
      dplyr::summarise(n = dplyr::n(), mean_shap = mean(shap_logodds),
                       sd_shap = stats::sd(shap_logodds), sum_shap = sum(shap_logodds),
                       .groups = "drop")
  }
}
importance    <- dplyr::bind_rows(importance)
level_summary <- dplyr::bind_rows(level_summary)
utils::write.csv(importance,    here::here("outputs", "shap_importance.csv"),    row.names = FALSE)
utils::write.csv(level_summary, here::here("outputs", "shap_level_summary.csv"), row.names = FALSE)
print(as.data.frame(importance), row.names = FALSE)

# =============================================================================
# 3. Plot functions
# =============================================================================
xlab_shap <- function(cls) paste0("SHAP value (log-odds, ", CLASS_LABEL[[cls]], ")")

## Importance: mean |SHAP| per feature
plot_importance <- function(sv_cls, cls, file) {
  S <- shapviz::get_shap_values(sv_cls)
  d <- tibble::tibble(feature = colnames(S), value = colMeans(abs(S)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = stats::reorder(feature, value))) +
    ggplot2::geom_col(fill = "#E69F00") +
    ggplot2::labs(x = paste0("Mean |SHAP value| (log-odds, ", CLASS_LABEL[[cls]], ")"), y = NULL) +
    theme_paper(24)
  ggplot2::ggsave(file, p, width = 12, height = 8)
}

## Dependence: one point per network; large point = level mean (not summed)
plot_dependence <- function(sv_cls, cls, feature, file) {
  d <- tibble::tibble(shap  = shapviz::get_shap_values(sv_cls)[, feature],
                      value = shapviz::get_feature_values(sv_cls)[[feature]])
  if (is.numeric(d$value)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = shap, y = value)) +
      ggplot2::geom_vline(xintercept = 0, colour = "red") +
      ggplot2::geom_point(alpha = 0.6, size = 2.5) +
      ggplot2::labs(x = xlab_shap(cls), y = feature)
  } else {
    d <- d %>% dplyr::add_count(value) %>%
      dplyr::mutate(label = sprintf("%s (n = %d)", tools::toTitleCase(as.character(value)), n))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = shap, y = stats::reorder(label, shap, FUN = mean))) +
      ggplot2::geom_vline(xintercept = 0, colour = "red") +
      ggplot2::geom_jitter(height = 0.2, width = 0, alpha = 0.4, size = 2) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 6) +
      ggplot2::labs(x = xlab_shap(cls), y = NULL, title = feature)
  }
  ggplot2::ggsave(file, p + theme_paper(24), width = 12, height = 8)
}

## Waterfall: one network
plot_waterfall <- function(sv_cls, cls, row, file) {
  p <- shapviz::sv_waterfall(sv_cls, row_id = row) +
    ggplot2::labs(x = paste0("Log-odds (", CLASS_LABEL[[cls]], ")")) + theme_paper(24)
  ggplot2::ggsave(file, p, width = 12, height = 8)
}

## Two-way: points at each pair of levels, coloured by the SHAP interaction value
plot_twoway <- function(sv_cls, x, y, file) {
  p <- shapviz::sv_dependence2D(sv_cls, x = x, y = y, interactions = TRUE) +
    ggplot2::labs(colour = "SHAP interaction value (log-odds)") + theme_paper(24) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom",
                   legend.justification = "center",
                   legend.box.just = "center",
                   legend.margin = ggplot2::margin(t = 5, r = 20, b = 5, l = 20)) +
    ## title above a wide, centred colour bar so the end labels are not cut off
    ggplot2::guides(colour = ggplot2::guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = ggplot2::unit(10, "cm"), barheight = ggplot2::unit(0.6, "cm")))
  ggplot2::ggsave(file, p, width = 12, height = 9)
}

# =============================================================================
# 4. Figures. File names are those used in the manuscript.
# =============================================================================
m2_geo <- shap_no_species$Spatial;  m2_sf <- shap_no_species$`Scale-Free`
m1_geo <- shap_species$Spatial;     m1_sf <- shap_species$`Scale-Free`

## Main text, Model 2, geometric graph: Figure 4 (importance), Figure 5 (dependence)
plot_importance(m2_geo, "Spatial", fig("metavarimp_sp-new.pdf"))
plot_dependence(m2_geo, "Spatial", "Data collection",      fig("data_collectn-new.dep2.pdf"))
plot_dependence(m2_geo, "Spatial", "Data duration (days)", fig("data_duratn-new.dep2.pdf"))
plot_dependence(m2_geo, "Spatial", "Captive",              fig("captive-new.dep2.pdf"))
plot_dependence(m2_geo, "Spatial", "Interaction type",     fig("interac_type-new.dep2.pdf"))

## Supplement, Model 2, geometric graph: waterfall
## (build without shuffling so WATERFALL_ROW refers to metadata_train2.rds)
m2_geo_unshuffled <- build_shap(here::here("models", "booster_portable_meta2.json"),
                                here::here("models", "preproc_trained_meta2.rds"),
                                here::here("models", "metadata_train2.rds"),
                                SHAP_NAMES_NO_SPECIES, shuffle = FALSE)$Spatial
plot_waterfall(m2_geo_unshuffled, "Spatial", WATERFALL_ROW, fig("sp.meta.waterfall.pdf"))

## Supplement, Model 1, geometric graph
plot_importance(m1_geo, "Spatial", fig("metavarimp_sp_species.pdf"))
plot_dependence(m1_geo, "Spatial", "Species",              fig("species.dep.pdf"))
plot_dependence(m1_geo, "Spatial", "Data collection",      fig("data_collectn.dep_species.pdf"))
plot_dependence(m1_geo, "Spatial", "Data duration (days)", fig("data_duratn.dep_species.pdf"))
plot_dependence(m1_geo, "Spatial", "Captive",              fig("captive.dep_species.pdf"))
plot_dependence(m1_geo, "Spatial", "Interaction type",     fig("interaction_type_species.pdf"))

## Supplement, Model 1, scale-free
plot_importance(m1_sf, "Scale-Free", fig("metavarimp_sf_species.pdf"))
plot_dependence(m1_sf, "Scale-Free", "Species",              fig("species.dep_sf.pdf"))
plot_dependence(m1_sf, "Scale-Free", "Data collection",      fig("data_collectn_SF_species.pdf"))
plot_dependence(m1_sf, "Scale-Free", "Data duration (days)", fig("data_duratn_sf_species.pdf"))
plot_dependence(m1_sf, "Scale-Free", "Captive",              fig("captive_sf_species.pdf"))
plot_dependence(m1_sf, "Scale-Free", "Interaction type",     fig("interctn_SF_species.pdf"))

## Supplement, Model 2, scale-free (duration panel has its own file)
plot_importance(m2_sf, "Scale-Free", fig("metavarimp_sf.pdf"))
plot_dependence(m2_sf, "Scale-Free", "Data collection",      fig("data_collectn_SF.pdf"))
plot_dependence(m2_sf, "Scale-Free", "Data duration (days)", fig("data_duratn_SF.pdf"))
plot_dependence(m2_sf, "Scale-Free", "Captive",              fig("cap_SF.pdf"))
plot_dependence(m2_sf, "Scale-Free", "Interaction type",     fig("interctn_SF.pdf"))

## Supplement, two-way plots, Model 2, geometric graph
plot_twoway(m2_geo, "Data duration (days)", "Interaction type", fig("dataduratn_interactn.pdf"))
plot_twoway(m2_geo, "Data duration (days)", "Data collection",  fig("dataduratn_collectn.pdf"))
plot_twoway(m2_geo, "Data duration (days)", "Captive",          fig("dataduratn_captive.pdf"))
plot_twoway(m2_geo, "Data collection",      "Interaction type", fig("datacollectn_interactn.pdf"))
plot_twoway(m2_geo, "Data collection",      "Captive",          fig("datacollectn_captive.pdf"))
plot_twoway(m2_geo, "Interaction type",     "Captive",          fig("Interaction_captive.pdf"))

## Supplement, two-way plots, Model 1, geometric graph
plot_twoway(m1_geo, "Species",          "Data duration (days)", fig("species_duration_interactn_sp.pdf"))
plot_twoway(m1_geo, "Species",          "Captive",              fig("species_Captive_interactn_sp.pdf"))
plot_twoway(m1_geo, "Species",          "Data collection",      fig("species_DataColl_interactn_sp.pdf"))
plot_twoway(m1_geo, "Species",          "Interaction type",     fig("species_intrc_interactn_sp.pdf"))
plot_twoway(m1_geo, "Interaction type", "Data collection",      fig("DataColl_interactn_sp.pdf"))

## Supplement, two-way plots, Model 2, scale-free
plot_twoway(m2_sf, "Data collection", "Interaction type",     fig("DataColl_interactn_sf.pdf"))
plot_twoway(m2_sf, "Data collection", "Class",                fig("DataColl_class_sf.pdf"))
plot_twoway(m2_sf, "Data collection", "Data duration (days)", fig("DataColl_duration_sf.pdf"))
plot_twoway(m2_sf, "Data collection", "Captive",              fig("DataColl_captive_sf.pdf"))

## Supplement, two-way plots, Model 1, scale-free
plot_twoway(m1_sf, "Species",          "Data collection",  fig("Species_collection_sf.pdf"))
plot_twoway(m1_sf, "Species",          "Interaction type", fig("Species_interaction_sf.pdf"))
plot_twoway(m1_sf, "Species",          "Data duration (days)", fig("Species_duration_sf.pdf"))
plot_twoway(m1_sf, "Data collection",  "Interaction type", fig("Datacollectn_interaction_species_sf.pdf"))
plot_twoway(m1_sf, "Data collection",  "Captive",          fig("Datacollectn_captive_species_sf.pdf"))
plot_twoway(m1_sf, "Interaction type", "Captive",          fig("Interaction_captive_species_sf.pdf"))

message("07_interpretation_shap_hstats.R complete.")