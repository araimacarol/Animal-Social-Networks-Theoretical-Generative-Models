# =============================================================================
# 05_group_level_community_detection.R
# -----------------------------------------------------------------------------
# Group-level (multilevel) analysis.
#
# We:
#   (a) Do not classify very small communities -> minimum group size of 8 nodes
#       (MIN_GROUP_SIZE in 00_setup.R). Field-defined badger social groups span
#       2-9 individuals, so a threshold of 8 keeps the largest genuine groups
#       and drops communities too small for stable structural features.
#   (b) To determne whether group-level results are stable across community-detection algorithms,
#       -> we run Infomap, leading-eigenvector, Louvain and Walktrap (plus the
#       prescribed badger grouping) and compare the predicted classes.
#
# SCOPE -- IMPORTANT
#   This analysis runs for badger, chimpanzee,and  giraffe, ant. 
#   It does not run on all empirical networks.
#
# Inputs : data/edges/<network>.edges, data/Badger_Groups.csv, recovered model
# Outputs: outputs/community_summary.csv, outputs/community_predictions.csv,
#          outputs/merged_community_predictions.csv,
#          outputs/community_method_stability.csv,
#          outputs/community_method_failures.csv
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

## Louvain and Infomap are stochastic; seed immediately before the run so the
## partitions reproduce.
set.seed(123)

model <- load_main_model()

## ---- The 41 multilevel networks ---------------------------------------------
MULTILEVEL_NETWORKS <- c(
  "Badger_Autumn1", "Badger_Autumn2", "Badger_Autumn3",
  "Badger_Autumn4", "Badger_Spring1", "Badger_Spring2",
  "Badger_Spring3", "Badger_Summer1", "Badger_Summer2",
  "Badger_Summer3", "Badger_Summer4", "Badger_Winter2",
  "Badger_Winter3", "Badger_Winter4", "Chimp_AprSRI_G_Matrix",
  "Chimp_Apr_Prob5_MatrixNoNAs", "Chimp_Aug_Prob5_MatrixNoNAs", "Chimp_Dec_Prob5_MatrixNoNAs",
  "Chimp_FebSRI_G_Matrix", "Chimp_Feb_Prob5_MatrixNoNAs", "Chimp_JanSRI_G_Matrix",
  "Chimp_Jan_Prob5_MatrixNoNAs", "Chimp_JulSRI_G_Matrix", "Chimp_Jul_Prob5_MatrixNoNAs",
  "Chimp_Jun_Prob5_MatrixNoNAs", "Chimp_MarSRI_G_Matrix", "Chimp_Mar_Prob5_MatrixNoNAs",
  "Chimp_MaySRI_G_Matrix", "Chimp_May_Prob5_MatrixNoNAs", "Chimp__AugSRI_G_Matrix",
  "Chimp__DecSRI_G_Matrix", "Chimp__JunSRI_G_Matrix", "Giraffe_edgelist",
  "Giraffe_group1", "Giraffe_group2", "insecta-ant-colony1",
  "insecta-ant-colony2", "insecta-ant-colony3", "insecta-ant-colony4",
  "insecta-ant-colony5", "insecta-ant-colony6"
)

## Resolve against the real file names in data/Empirical-Network-Edges (see
## resolve_edge_files() in 01_functions.R -- handles the ".csv .edges" suffix).
edge_files <- resolve_edge_files(MULTILEVEL_NETWORKS, dir = dir_edges)
stopifnot(length(edge_files) == length(MULTILEVEL_NETWORKS))

message("Running community detection on ", length(edge_files), " multilevel networks.")

## ---- Group sizes + predicted class, per network, per method -----------------
cd <- summarize_group_sizes_and_predictions(
  edge_files     = edge_files,
  model          = model,
  methods        = c("walktrap", "louvain", "leading_eigen", "infomap", "prescribed"),
  min_group_size = MIN_GROUP_SIZE
)

utils::write.csv(cd$size_summary,
                 here::here("outputs", "community_summary.csv"), row.names = FALSE)
utils::write.csv(cd$predictions,
                 here::here("outputs", "community_predictions.csv"), row.names = FALSE)

## ---- Algorithm failures (ARPACK non-convergence etc.) -----------------------
## Reported, not silently dropped: leading-eigenvector can fail to converge, and
## the partition igraph returns in that case is not a valid result.
if (nrow(cd$failures) > 0) {
  cat("\n### COMMUNITY-DETECTION FAILURES (excluded from all results) ###\n")
  print(cd$failures %>% dplyr::count(method, reason, name = "n_networks"))
  utils::write.csv(cd$failures,
                   here::here("outputs", "community_method_failures.csv"), row.names = FALSE)
} else {
  message("All methods converged on all networks.")
}

## ---- Merge into the single table --------------------------
merged_table <- merge_community_predictions(cd)
utils::write.csv(merged_table,
                 here::here("outputs", "merged_community_predictions.csv"), row.names = FALSE)

## ---- Cross-method stability summary (reported in the Results) ---------------
stability <- merged_table %>%
  dplyr::filter(method != "prescribed", !is.na(predicted_class_label)) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    n_networks    = dplyr::n_distinct(graph_name),
    n_communities = dplyr::n(),
    n_geometric   = sum(predicted_class_label == "Geometric Graph"),
    pct_geometric = round(100 * n_geometric / n_communities, 1),
    .groups = "drop")

cat("\n### CROSS-METHOD STABILITY (communities >= ", MIN_GROUP_SIZE, " nodes) ###\n", sep = "")
print(stability)
utils::write.csv(stability,
                 here::here("outputs", "community_method_stability.csv"), row.names = FALSE)

## ---- Rows per method: sanity check against the published run ----------------
cat("\n### ROWS PER METHOD (published run: walktrap 62, louvain 76,",
    "leading_eigen 76, infomap 54, prescribed 21) ###\n")
print(merged_table %>% dplyr::count(method, name = "n_rows"))
cat("Total rows:", nrow(merged_table), " (published: 289)\n")

message("05_group_level_community_detection.R complete.")