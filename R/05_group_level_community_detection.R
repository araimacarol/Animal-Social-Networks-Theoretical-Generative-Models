# =============================================================================
# 05_group_level_community_detection.R
# -----------------------------------------------------------------------------
# Group-level (multilevel) analysis.
#
# We:
#   (a) Do not classify very small groups -> minimum group size of 8 nodes
#       (MIN_GROUP_SIZE in 00_setup.R). Field-defined badger social groups span
#       2-9 individuals, so a threshold of 8 keeps the largest genuine groups
#       and drops groups too small for stable structural features.
#   (b) Test whether group-level results are stable across community-detection
#       algorithms -> run Walktrap, Louvain, leading-eigenvector and Infomap and
#       compare the predicted classes.
#
# SCOPE -- IMPORTANT
#   * Badger networks are EXCLUDED from the algorithm comparison: their groups
#     are field-defined (Silk, Weber et al. 2018), not detected, so running
#     detection algorithms on them cannot test detection stability. They are
#     handled separately using their prescribed groups.
#   * The ant set includes the six aggregate colony networks AND the two daily
#     colony-1 networks with the highest modularity (day20, day05). The
#     aggregate colony networks have negligible modularity (< 0.01), so any
#     partition of them is an arbitrary division of an essentially unstructured
#     graph; the daily networks carry real community structure and are the
#     meaningful test for ants. Selection is by MODULARITY, not by predicted
#     class -- selecting on the outcome would be circular.
#
# Inputs : data/Empirical-Network-Edges/<network>.edges, data/Badger_Groups.csv,
#          recovered first-stage model
# Outputs: outputs/community_summary.csv, outputs/community_predictions.csv,
#          outputs/merged_community_predictions.csv,
#          outputs/community_method_stability.csv,
#          outputs/community_network_modularity.csv,
#          outputs/community_method_failures.csv,
#          outputs/badger_prescribed_groups.csv,
#          outputs/table_s6_representative.csv
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

## Louvain and Infomap are stochastic; seed immediately before the run so the
## partitions reproduce.
set.seed(123)

model <- load_main_model()

## ---- Networks entering the ALGORITHM COMPARISON (detected groups) -----------
CHIMP_NETWORKS <- c(
  "Chimp_AprSRI_G_Matrix", "Chimp_Apr_Prob5_MatrixNoNAs",
  "Chimp_Aug_Prob5_MatrixNoNAs", "Chimp_Dec_Prob5_MatrixNoNAs",
  "Chimp_FebSRI_G_Matrix", "Chimp_Feb_Prob5_MatrixNoNAs",
  "Chimp_JanSRI_G_Matrix", "Chimp_Jan_Prob5_MatrixNoNAs",
  "Chimp_JulSRI_G_Matrix", "Chimp_Jul_Prob5_MatrixNoNAs",
  "Chimp_Jun_Prob5_MatrixNoNAs", "Chimp_MarSRI_G_Matrix",
  "Chimp_Mar_Prob5_MatrixNoNAs", "Chimp_MaySRI_G_Matrix",
  "Chimp_May_Prob5_MatrixNoNAs", "Chimp__AugSRI_G_Matrix",
  "Chimp__DecSRI_G_Matrix", "Chimp__JunSRI_G_Matrix"
)

GIRAFFE_NETWORKS <- c("Giraffe_edgelist", "Giraffe_group1", "Giraffe_group2")

## Aggregate colonies (near-zero modularity) + the two highest-modularity daily
## colony-1 networks. See SCOPE above for why both are included.
ANT_NETWORKS <- c(
  paste0("insecta-ant-colony", 1:6),
  "insecta-ant-colony1-day20",   # modularity 0.1604
  "insecta-ant-colony1-day05"    # modularity 0.1399
)

DETECTION_NETWORKS <- c(CHIMP_NETWORKS, GIRAFFE_NETWORKS, ANT_NETWORKS)

## ---- Badger networks: PRESCRIBED (field-defined) groups only ----------------
BADGER_NETWORKS <- c(
  "Badger_Autumn1", "Badger_Autumn2", "Badger_Autumn3", "Badger_Autumn4",
  "Badger_Spring1", "Badger_Spring2", "Badger_Spring3",
  "Badger_Summer1", "Badger_Summer2", "Badger_Summer3", "Badger_Summer4",
  "Badger_Winter2", "Badger_Winter3", "Badger_Winter4"
)

detection_files <- resolve_edge_files(DETECTION_NETWORKS, dir = dir_edges)
badger_files    <- resolve_edge_files(BADGER_NETWORKS,    dir = dir_edges)
stopifnot(length(detection_files) == length(DETECTION_NETWORKS),
          length(badger_files)    == length(BADGER_NETWORKS))

message("Algorithm comparison: ", length(detection_files),
        " networks (", length(CHIMP_NETWORKS), " chimpanzee, ",
        length(GIRAFFE_NETWORKS), " giraffe, ", length(ANT_NETWORKS), " ant).")
message("Prescribed groups:    ", length(badger_files), " badger networks.")

## ---- Per-network modularity (context for the stability figures) -------------
## Walktrap modularity of the largest component. Networks with modularity near
## zero have no community structure for any algorithm to recover, so their
## "groups" are arbitrary. Reported so this can be stated rather than assumed.
network_modularity <- purrr::map_dfr(detection_files, function(f) {
  gname <- tools::file_path_sans_ext(basename(f))
  g <- tryCatch(largest_component(read_edges_graph(f)), error = function(e) NULL)
  if (is.null(g) || igraph::vcount(g) < MIN_GROUP_SIZE) return(NULL)
  tibble::tibble(
    graph_name = gname,
    n_nodes    = as.integer(igraph::vcount(g)),
    modularity = round(igraph::modularity(g,
                    igraph::membership(igraph::cluster_walktrap(g))), 4))
})

utils::write.csv(network_modularity,
                 here::here("outputs", "community_network_modularity.csv"), row.names = FALSE)

cat("\n### NETWORKS WITH NEGLIGIBLE MODULARITY (< 0.01) ###\n")
print(network_modularity %>% dplyr::filter(modularity < 0.01))

## ---- 1. Algorithm comparison (detected groups, no badger) -------------------
cd_detect <- summarize_group_sizes_and_predictions(
  edge_files     = detection_files,
  model          = model,
  methods        = c("walktrap", "louvain", "leading_eigen", "infomap"),
  min_group_size = MIN_GROUP_SIZE
)

## ---- 2. Badger: prescribed groups only --------------------------------------
cd_badger <- summarize_group_sizes_and_predictions(
  edge_files     = badger_files,
  model          = model,
  methods        = "prescribed",
  min_group_size = MIN_GROUP_SIZE
)

utils::write.csv(cd_detect$size_summary,
                 here::here("outputs", "community_summary.csv"), row.names = FALSE)
utils::write.csv(cd_detect$predictions,
                 here::here("outputs", "community_predictions.csv"), row.names = FALSE)

badger_table <- merge_community_predictions(cd_badger)
utils::write.csv(badger_table,
                 here::here("outputs", "badger_prescribed_groups.csv"), row.names = FALSE)

## ---- Algorithm failures (ARPACK non-convergence etc.) -----------------------
all_failures <- dplyr::bind_rows(cd_detect$failures, cd_badger$failures)
if (nrow(all_failures) > 0) {
  cat("\n### COMMUNITY-DETECTION FAILURES (excluded from all results) ###\n")
  print(all_failures %>% dplyr::count(method, reason, name = "n_networks"))
  utils::write.csv(all_failures,
                   here::here("outputs", "community_method_failures.csv"), row.names = FALSE)
} else {
  message("All methods converged on all networks.")
}

## ---- Merged table: detected groups + badger prescribed groups ---------------
detect_table <- merge_community_predictions(cd_detect)
merged_table <- dplyr::bind_rows(detect_table, badger_table)
utils::write.csv(merged_table,
                 here::here("outputs", "merged_community_predictions.csv"), row.names = FALSE)

## ---- Cross-method stability (detected groups only) --------------------------
stability_all <- detect_table %>%
  dplyr::filter(!is.na(predicted_class_label)) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    n_networks    = dplyr::n_distinct(graph_name),
    n_groups      = dplyr::n(),
    n_geometric   = sum(predicted_class_label == "Geometric Graph"),
    pct_geometric = round(100 * n_geometric / n_groups, 1),
    .groups = "drop") %>%
  dplyr::mutate(subset = "all networks")

## Same figures restricted to networks with real community structure. Reported
## alongside, not instead of, the full set.
structured <- network_modularity$graph_name[network_modularity$modularity >= 0.01]

stability_structured <- detect_table %>%
  dplyr::filter(!is.na(predicted_class_label),
                tools::file_path_sans_ext(graph_name) %in% structured |
                  graph_name %in% structured) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    n_networks    = dplyr::n_distinct(graph_name),
    n_groups      = dplyr::n(),
    n_geometric   = sum(predicted_class_label == "Geometric Graph"),
    pct_geometric = round(100 * n_geometric / n_groups, 1),
    .groups = "drop") %>%
  dplyr::mutate(subset = "modularity >= 0.01")

stability <- dplyr::bind_rows(stability_all, stability_structured)

cat("\n### CROSS-METHOD STABILITY (groups >= ", MIN_GROUP_SIZE, " nodes) ###\n", sep = "")
print(stability)
utils::write.csv(stability,
                 here::here("outputs", "community_method_stability.csv"), row.names = FALSE)

## ---- Badger summary (reported separately) -----------------------------------
badger_summary <- badger_table %>%
  dplyr::filter(!is.na(predicted_class_label)) %>%
  dplyr::summarise(
    n_networks    = dplyr::n_distinct(graph_name),
    n_groups      = dplyr::n(),
    n_geometric   = sum(predicted_class_label == "Geometric Graph"),
    pct_geometric = round(100 * n_geometric / n_groups, 1))

cat("\n### BADGER, FIELD-DEFINED GROUPS (groups >= ", MIN_GROUP_SIZE, " nodes) ###\n", sep = "")
print(badger_summary)

## ---- Representative subset for Table S6 -------------------------------------
## One chimpanzee, one giraffe, and the ant aggregate plus its two daily
## networks, so the ant rows show both the unstructured aggregate and the
## structured daily networks.
S6_NETWORKS <- c("Chimp_Dec_Prob5_MatrixNoNAs", "Giraffe_group1",
                 "insecta-ant-colony1",
                 "insecta-ant-colony1-day05", "insecta-ant-colony1-day20")

table_s6 <- detect_table %>%
  dplyr::mutate(base_name = trimws(tools::file_path_sans_ext(graph_name))) %>%
  dplyr::filter(base_name %in% S6_NETWORKS, !is.na(predicted_class_label)) %>%
  dplyr::left_join(network_modularity, by = c("base_name" = "graph_name")) %>%
  dplyr::select(network = base_name, method, total_nodes_in_network, modularity,
                group_id, n_nodes_from_passed_group, predicted_class_label) %>%
  dplyr::arrange(match(network, S6_NETWORKS), method, group_id)

cat("\n### TABLE S6 REPRESENTATIVE SUBSET ###\n")
print(as.data.frame(table_s6), row.names = FALSE)
utils::write.csv(table_s6,
                 here::here("outputs", "table_s6_representative.csv"), row.names = FALSE)

## ---- Rows per method: scope check -------------------------------------------
cat("\n### ROWS PER METHOD ###\n")
print(merged_table %>% dplyr::count(method, name = "n_rows"))
cat("Total rows:", nrow(merged_table), "\n")

message("05_group_level_community_detection.R complete.")
