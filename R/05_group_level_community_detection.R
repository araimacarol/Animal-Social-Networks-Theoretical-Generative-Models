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
#   (b) To determine whether group-level results are stable across community-
#       detection algorithms, we compare leading-eigenvector, Louvain and
#       Walktrap.
#       Infomap returned a single community spanning the whole network for every
#       chimpanzee, giraffe and ant network, so it does not subdivide these
#       networks and is not used in the comparison. It is
#       still RUN, unchanged, because Louvain is stochastic and draws from the
#       same random-number stream: removing Infomap from the run would change
#       every Louvain partition. Its rows are dropped from all outputs.
#   (c) The cross-method comparison uses the chimpanzee, giraffe and ant
#       networks (badger groups are already field-defined), overall and restricted to
#       networks with detectable community structure (modularity >= 0.01;
#       walktrap modularity of the whole network).
#
# SCOPE -- IMPORTANT
#   This analysis runs for badger, chimpanzee, giraffe and ant networks
#   (six aggregate colony networks plus the colony-1 day-5 and day-20 networks).
#   It does not run on all empirical networks.
#
# Inputs : data/Empirical-Network-Edges/ <network>.edges, data/Badger_Groups.csv, recovered model
# Outputs: outputs/community_summary.csv, outputs/community_predictions.csv,
#          outputs/merged_community_predictions.csv,
#          outputs/community_method_stability.csv,
#          outputs/community_network_modularity.csv,
#          outputs/community_method_failures.csv
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

## Louvain and Infomap are stochastic and share R's random-number stream; seed
## immediately before the run so the partitions reproduce.
set.seed(123)

model <- load_main_model()

## ---- The 43 multilevel networks (41 + 2 ant daily networks) ---------------------------------------------
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
  "insecta-ant-colony5", "insecta-ant-colony6",
  ## Two colony-1 daily networks with the highest modularity (selected on
  ## modularity, not on predicted class). Appended LAST on purpose: Louvain and
  ## Infomap draw from the random-number stream network by network, so adding
  ## them at the end leaves every partition of the 41 networks above unchanged.
  "insecta-ant-colony1-day05", "insecta-ant-colony1-day20"
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

## Drop Infomap rows (see header) and clean network names: file names carry a
## trailing ".csv " (e.g. "Chimp_AprSRI_G_Matrix.csv "), which would break joins.
clean_name <- function(x) sub("\\.csv$", "", trimws(x))
cd$size_summary <- cd$size_summary %>%
  dplyr::filter(method != "infomap") %>% dplyr::mutate(graph_name = clean_name(graph_name))
cd$predictions  <- cd$predictions %>%
  dplyr::filter(method != "infomap") %>% dplyr::mutate(graph_name = clean_name(graph_name))
utils::write.csv(cd$size_summary,
                 here::here("outputs", "community_summary.csv"), row.names = FALSE)
utils::write.csv(cd$predictions,
                 here::here("outputs", "community_predictions.csv"), row.names = FALSE)
## (failures, if any, are reported below for all methods actually run)

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

## ---- Merge into the single table (Infomap rows dropped; see header) ---------
merged_table <- merge_community_predictions(cd)
utils::write.csv(merged_table,
                 here::here("outputs", "merged_community_predictions.csv"), row.names = FALSE)

## ---- Network modularity (walktrap, whole network; deterministic) -------------
modularity_tab <- dplyr::bind_rows(lapply(edge_files, function(f) {
  ## same graph construction as summarize_group_sizes_and_predictions()
  df <- utils::read.table(f, sep = " ", fill = TRUE, header = FALSE)
  g  <- igraph::graph_from_data_frame(as.matrix(df), directed = FALSE)
  g  <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g  <- largest_component(igraph::delete_vertices(g, which(igraph::degree(g) == 0)))
  tibble::tibble(graph_name = clean_name(tools::file_path_sans_ext(basename(f))),
                 modularity = igraph::modularity(igraph::cluster_walktrap(g)))
}))
stopifnot(all(unique(merged_table$graph_name) %in% modularity_tab$graph_name))
utils::write.csv(modularity_tab, here::here("outputs", "community_network_modularity.csv"),
                 row.names = FALSE)

## ---- Cross-method stability summary (reported in the Results) ---------------
## Comparison set: chimpanzee, giraffe and ant networks (badgers excluded).
stab_fun <- function(d) d %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    n_networks       = dplyr::n_distinct(graph_name),
    n_networks_split = dplyr::n_distinct(graph_name[n_groups_detected >= 2]),
    n_communities    = sum(!is.na(predicted_class_label)),   # classified (>= MIN_GROUP_SIZE)
    n_geometric      = sum(predicted_class_label == "Geometric Graph", na.rm = TRUE),
    pct_geometric    = round(100 * n_geometric / n_communities, 1),
    .groups = "drop")
comp <- merged_table %>%
  dplyr::filter(method != "prescribed", !grepl("^Badger", graph_name)) %>%
  dplyr::left_join(modularity_tab, by = "graph_name")
stability <- dplyr::bind_rows(
  dplyr::mutate(stab_fun(comp), subset = "all comparison networks", .before = 1),
  dplyr::mutate(stab_fun(dplyr::filter(comp, modularity >= 0.01)),
                subset = "modularity >= 0.01", .before = 1))

cat("\n### CROSS-METHOD STABILITY (communities >= ", MIN_GROUP_SIZE, " nodes) ###\n", sep = "")
print(stability)
utils::write.csv(stability,
                 here::here("outputs", "community_method_stability.csv"), row.names = FALSE)

## ---- Rows per method ----------------------------------------------------------
print(merged_table %>% dplyr::count(method, name = "n_rows"))

message("05_group_level_community_detection.R complete.")