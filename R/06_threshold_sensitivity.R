# # =============================================================================
# # 06_threshold_sensitivity.R
# # -----------------------------------------------------------------------------
# # Binarizing the networks and taking the largest connected
# # component may make them look artificially "geometric". Does the predicted
# # class survive removing weak ties?
# #
# #The empirical network files are unweighted, so we cannot threshold true
# # interaction weights. Instead we derive a structural connection-strength proxy
# # (neighbourhood overlap / Jaccard of shared neighbours; optionally
# # 1/effective-resistance) and progressively delete the structurally weakest
# # edges, re-extracting the largest connected component and re-classifying at
# # each percentile of edge removal.
# #
# 
# #
# # Inputs : data/Networks-Threshold-Test/*.edges, recovered main model
# # Outputs: outputs/threshold_sensitivity_{raw,by_graph,wide}.csv
# # =============================================================================
# 
# source(here::here("R", "00_setup.R"))
# source(here::here("R", "01_functions.R"))
# 
# set.seed(123)
# 
# out <- threshold_sensitivity_test(
#   edges_dir     = here::here("data", "Networks-Threshold-Test"),
#   model         = load_main_model(),
#   weight_method = "adjacency",   # or "laplacian"
#   n_thresholds  = 10,
#   out_prefix    = here::here("outputs", "threshold_sensitivity")
# )
# 
# ## Per-graph summary: is the predicted class stable under pruning?
# print(out$by_graph, n = Inf)
# 
# ## Which networks changed class under any threshold?
# message("\nNetworks whose predicted class changed under pruning:")
# print(out$by_graph %>% dplyr::filter(!stable))
# 
# ## ---- Headline summary by predicted class (the main-text table) --------------
# summary_by_class <- out$results %>%
#   dplyr::filter(component_mode == "largest_component") %>%
#   dplyr::group_by(graph_name, baseline_class) %>%
#   dplyr::summarise(
#     inv_30 = all(matches_baseline[threshold_quant %in% c("0%", "10%", "20%", "30%")]),
#     inv_50 = all(matches_baseline[threshold_quant %in% c("0%", "10%", "20%", "30%", "40%", "50%")]),
#     .groups = "drop") %>%
#   dplyr::group_by(baseline_class) %>%
#   dplyr::summarise(
#     N = dplyr::n(),
#     pct_invariant_30 = round(100 * mean(inv_30, na.rm = TRUE), 1),
#     pct_invariant_50 = round(100 * mean(inv_50, na.rm = TRUE), 1),
#     .groups = "drop")
# 
# print(summary_by_class)
# utils::write.csv(summary_by_class,
#                  here::here("outputs", "threshold_sensitivity_by_class.csv"), row.names = FALSE)
# 
# message("06_threshold_sensitivity.R complete.")

# =============================================================================
# 06_threshold_sensitivity.R
# -----------------------------------------------------------------------------
# Binarising the networks and taking the largest connected
# component may make them look artificially "geometric". Does the predicted
# class survive removing weak ties?
#
# The empirical .edges files are unweighted, so we cannot threshold true
# interaction weights. Instead we derive a structural connection-strength proxy
# (neighbourhood overlap / Jaccard of shared neighbours; optionally
# 1/effective-resistance) and progressively delete the structurally weakest
# edges, re-extracting the largest connected component and re-classifying at
# each percentile of edge removal.
#
#
# Inputs : data/Networks-Threshold-Test/*.edges, recovered main model
# Outputs: outputs/threshold_sensitivity_{raw,by_graph,wide}.csv
# =============================================================================
source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))
set.seed(123)

out <- threshold_sensitivity_test(
  edges_dir     = here::here("data", "Networks-Threshold-Test"),
  model         = load_main_model(),
  weight_method = "adjacency",   # or "laplacian"
  n_thresholds  = 10,
  out_prefix    = here::here("outputs", "threshold_sensitivity")
)

## Per-graph summary: is the predicted class stable under pruning?
print(out$by_graph, n = Inf)

## Which networks changed class under any threshold?
message("\nNetworks whose predicted class changed under pruning:")
print(out$by_graph %>% dplyr::filter(!stable))

## ---- Headline summary by predicted class (the main-text table) --------------
## Invariance is assessed only over states in which the network remains
## classifiable; fragmentation is NOT counted as a class change. This matches
## the convention used by out$by_graph$stable and the wide file.
summary_by_class <- out$results %>%
  dplyr::filter(component_mode == "largest_component") %>%
  dplyr::group_by(graph_name, baseline_class) %>%
  dplyr::summarise(
    inv_30 = all(matches_baseline[threshold_quant %in% c("0%", "10%", "20%", "30%") &
                                    !is.na(predicted_class)]),
    inv_50 = all(matches_baseline[threshold_quant %in% c("0%", "10%", "20%", "30%", "40%", "50%") &
                                    !is.na(predicted_class)]),
    .groups = "drop") %>%
  dplyr::group_by(baseline_class) %>%
  dplyr::summarise(
    N = dplyr::n(),
    pct_invariant_30 = round(100 * mean(inv_30, na.rm = TRUE), 1),
    pct_invariant_50 = round(100 * mean(inv_50, na.rm = TRUE), 1),
    .groups = "drop")

print(summary_by_class)
utils::write.csv(summary_by_class,
                 here::here("outputs", "threshold_sensitivity_by_class.csv"), row.names = FALSE)

message("06_threshold_sensitivity.R complete.")
