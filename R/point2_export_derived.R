# =============================================================================
# point2_export_derived.R  --  RUN ONCE, LOCALLY (helper for PCA_Applicability.R)
# -----------------------------------------------------------------------------
# Reduces three large local files to the small derived files that
# PCA_Applicability.R needs for parts (a), (c) and (d). After this has run,
# the large files are no longer needed by the project.
#
# Before running: copy these three files into models/ (temporarily):
#   df_train_new.rds        309,499 pre-SMOTE training networks
#   final.fitted.model.rds  fitted workflow (its mold = the Figure 2 training rows)
#   final.emp.data.rds      636 empirical networks (Figure 2 feature set)
#
#   source(here::here("R", "00_setup.R"))
#   source(here::here("R", "point2_export_derived.R"))
#
# Writes to models/ (all small):
#   point2_training_feature_ranges.csv   min / max of each feature, training networks   (c)
#   point2_training_size_by_class.csv    size range per generative class               (a)
#   point2_fig2_pca.rds                  Figure 2 PCA: coordinates, variance, loadings (d)
# Keep final.emp.data.rds in models/ (small). Afterwards, delete df_train_new.rds
# and final.fitted.model.rds from models/.
# =============================================================================

M    <- here::here("models")
need <- file.path(M, c("df_train_new.rds", "final.fitted.model.rds", "final.emp.data.rds"))
if (!all(file.exists(need)))
  stop("Copy into models/ first:\n  ", paste(basename(need[!file.exists(need)]), collapse = "\n  "))

emp_fig <- readRDS(file.path(M, "final.emp.data.rds"))

## ---- (a), (c): training ranges from the pre-SMOTE training networks ---------
train_df   <- readRDS(file.path(M, "df_train_new.rds"))
feat_names <- intersect(setdiff(names(emp_fig), "graph_name"), names(train_df))

ranges <- tibble::tibble(feature   = feat_names,
                         train_min = sapply(train_df[feat_names], min, na.rm = TRUE),
                         train_max = sapply(train_df[feat_names], max, na.rm = TRUE))

cls_col  <- intersect(c("graph_name", "target", "GraphName"), names(train_df))[1]
size_tab <- train_df %>%
  dplyr::group_by(class = as.character(.data[[cls_col]])) %>%
  dplyr::summarise(n_networks = dplyr::n(), min_order = min(order), max_order = max(order),
                   n_order_lt_50 = sum(order < 50), .groups = "drop") %>%
  dplyr::bind_rows(dplyr::summarise(train_df, class = "All", n_networks = dplyr::n(),
                                    min_order = min(order), max_order = max(order),
                                    n_order_lt_50 = sum(order < 50)))

utils::write.csv(ranges,   file.path(M, "point2_training_feature_ranges.csv"), row.names = FALSE)
utils::write.csv(size_tab, file.path(M, "point2_training_size_by_class.csv"),  row.names = FALSE)
cat("Training networks:", nrow(train_df), "| features:", length(feat_names), "\n")
print(size_tab)
rm(train_df)

## ---- (d): the Figure 2 PCA (workflow mold + empirical, as published) --------
fit  <- readRDS(file.path(M, "final.fitted.model.rds"))
wf   <- tryCatch(tune::extract_workflow(fit), error = function(e) fit$.workflow[[1]])
mold <- hardhat::extract_mold(wf)$predictors %>% tibble::as_tibble()
rm(fit, wf)

combined <- dplyr::bind_rows(
  dplyr::mutate(mold, src = "sim"),
  emp_fig %>% dplyr::select(dplyr::all_of(colnames(mold))) %>% dplyr::mutate(src = "emp"))
n_sim <- nrow(mold); rm(mold)

pp <- recipes::recipe(~ ., data = combined) %>%
  recipes::update_role(src, new_role = "id") %>%
  recipes::step_normalize(recipes::all_predictors()) %>%
  recipes::step_pca(recipes::all_predictors(), num_comp = 2) %>%
  recipes::prep()
pc <- recipes::bake(pp, new_data = NULL)

v <- recipes::tidy(pp, number = 2, type = "variance") %>%
  dplyr::filter(terms == "percent variance", component <= 2) %>% dplyr::pull(value)

pca_out <- list(
  variance = c(PC1 = v[1], PC2 = v[2]),
  sim      = pc[pc$src == "sim", c("PC1", "PC2")],
  emp      = tibble::tibble(graph_name = as.character(emp_fig$graph_name), order = emp_fig$order,
                            PC1 = pc$PC1[pc$src == "emp"], PC2 = pc$PC2[pc$src == "emp"]),
  loadings = recipes::tidy(pp, number = 2, type = "coef") %>%
    dplyr::filter(component %in% c("PC1", "PC2")) %>%
    dplyr::select(terms, value, component))
saveRDS(pca_out, file.path(M, "point2_fig2_pca.rds"), compress = "xz")

## ---- Check against the published values -------------------------------------
hull <- pca_out$sim[grDevices::chull(pca_out$sim$PC1, pca_out$sim$PC2), ]
cov  <- mean(sp::point.in.polygon(pca_out$emp$PC1, pca_out$emp$PC2, hull$PC1, hull$PC2) > 0)
cat(sprintf("\nSimulated rows: %d (expect 474,235) | empirical: %d (expect 636)\n", n_sim, nrow(pca_out$emp)))
cat(sprintf("PC1 = %.1f%%, PC2 = %.1f%% (expect 38.3 / 22.3) | convex hull = %.1f%% (expect 95.9)\n",
            v[1], v[2], 100 * cov))
cat("\nWritten to models/. You can now delete df_train_new.rds and final.fitted.model.rds from models/.\n")
