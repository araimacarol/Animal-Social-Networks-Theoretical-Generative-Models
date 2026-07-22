# =============================================================================
# 02_first_stage_and_metadata.R
# -----------------------------------------------------------------------------
# 1. Classify every empirical network with the RECOVERED first-stage model.
# 2. Join predicted classes to the network metadata.
# 3. Write the analysis CSVs and per-category sample-size tables.
#
# MODEL ROUTE -- IMPORTANT
#   We do NOT use tune::extract_workflow(xgb_final_trained_model.rds) + predict().
#   That workflow carries the xgb.Booster that throws
#     "object is corrupted or is from an incompatible XGBoost version"
#   under current xgboost. We use the recovered version-portable JSON booster plus the
#   trained recipe (load_main_model() / predict_main_model() in 01_functions.R),
#   which is the same route scripts 05 and 06 use. Numerically identical to the
#   original workflow predict(): same booster, same trained recipe.
#
# Inputs : data/GraphFeatOnAllAnimNets.csv
#          data/Network_repository_metaData_combined.csv
#          models/booster_portable_mainmodel.json
#          models/preproc_trained_mainmodel.rds
# Outputs: outputs/Predicted-Class-of-Animal-Social-network.csv
#          outputs/meta_combined_df.csv, outputs/sample_sizes_*.csv
#          figures/classification_by_data_collection.png
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

set.seed(9456)

## ---- 1. First-stage predictions on empirical networks -----------------------

model <- load_main_model()   # portable booster + trained recipe + CLASS_ORDER

emp_data <- utils::read.csv(here::here("data", "GraphFeatOnAllAnimNets.csv")) %>%
  dplyr::rename(graph_name = GraphNames) %>%
  dplyr::select(graph_name, order, edges, mean_eccentr, mean_path_length, graph_energy,
                modularity, diameter, betw_centr, transitivity, spectral_radius,
                eigen_centr, deg_centr, mean_degree, minCut, FiedlerValue,
                Normalized_FiedlerValue, closeness_centr, deg_assort_coef) %>%
  janitor::clean_names() %>%
  dplyr::filter(order > 10) %>%
  # read.csv yields integers for order/edges/diameter/min_cut; step_impute_median()
  # in current recipes errors on integer columns.
  dplyr::mutate(dplyr::across(dplyr::where(is.integer), as.double)) %>%
  tibble::as_tibble()

message("Empirical networks after order > 10 filter: ", nrow(emp_data))
print(colSums(is.na(emp_data)))

# Empirical median imputation, exactly as in the original script (monolith L363-366).
# Medians are taken from the empirical data, not the training set -- keep it this
# way or the published first-stage predictions will not reproduce.
emp_prepped <- recipes::recipe(graph_name ~ ., data = emp_data) %>%
  recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
  recipes::prep() %>%
  recipes::juice()

stopifnot(!anyNA(emp_prepped))

# predict_main_model() bakes the TRAINED recipe, drops graph_name, and scores
# with the native booster. Row order is preserved.
preds <- predict_main_model(model, emp_prepped)
stopifnot(nrow(preds) == nrow(emp_prepped))

df_emp <- tibble::tibble(
  target            = as.character(emp_prepped$graph_name),
  Predicted_classes = factor(preds$predicted_class, levels = CLASS_ORDER)
)

utils::write.csv(df_emp,
                 here::here("outputs", "Predicted-Class-of-Animal-Social-network.csv"),
                 row.names = FALSE)

print(table(df_emp$Predicted_classes))

## ---- 2. Join to metadata ----------------------------------------------------

metadata <- utils::read.csv(here::here("data", "Network_repository_metaData_combined.csv")) %>%
  dplyr::select(-dplyr::any_of(c("Source", "Citation", "Edge_weight",
                                 "Time_resolution_secs"))) %>%
  dplyr::mutate(Graph.Name = gsub("_", "-", Graph.Name))

df_emp$target <- gsub("_", "-", df_emp$target)

meta_ordered <- metadata[match(df_emp$target, metadata$Graph.Name), ]
message("Networks with no metadata match: ", sum(is.na(meta_ordered$Graph.Name)))

combined_df <- dplyr::bind_cols(meta_ordered, df_emp) %>%
  dplyr::mutate(
    Predicted_classes = unname(CLASS_DISPLAY[as.character(Predicted_classes)])
  )

utils::write.csv(combined_df, here::here("outputs", "meta_combined_df.csv"),
                 row.names = FALSE)

## ---- 3. Per-category sample-size tables ----------

sample_sizes <- list(
  by_class       = combined_df %>% dplyr::count(Class, name = "n_networks"),
  by_interaction = combined_df %>% dplyr::count(Class, Interaction_type, name = "n_networks") %>%
    dplyr::arrange(Class, dplyr::desc(n_networks)),
  by_collection  = combined_df %>% dplyr::count(Class, Data_collection, name = "n_networks") %>%
    dplyr::arrange(Class, dplyr::desc(n_networks)),
  by_captive     = combined_df %>% dplyr::count(Class, Captive, name = "n_networks") %>%
    dplyr::arrange(Class, dplyr::desc(n_networks))
)
for (nm in names(sample_sizes)) {
  utils::write.csv(sample_sizes[[nm]],
                   here::here("outputs", paste0("sample_sizes_", nm, ".csv")),
                   row.names = FALSE)
}

## ---- 4. Main stacked-bar figure ---------------------------------------------

p_taxa <- metadata_plot_taxa(combined_df, facet_column = "Class",
                             group_by_column = "Data_collection")
ggplot2::ggsave(here::here("figures", "classification_by_data_collection.png"),
                p_taxa, width = 14, height = 6)

message("02_first_stage_and_metadata.R complete.")