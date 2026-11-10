# =============================================================================
# 06_threshold_sensitivity.R
# -----------------------------------------------------------------------------
# Sensitivity of the predicted generative-model class to network topology.
#
# All networks were binarised (no edge weights are used anywhere in this
# study). Thresholding observed interaction weights is therefore not possible.
# Instead, each edge is scored with a structural connection-strength proxy, the
# neighbourhood overlap (Jaccard index of shared neighbours), and the
# structurally weakest edges are removed in steps of 10 percentiles; the
# network is re-classified after each step. This is a TOPOLOGICAL perturbation
# test. It does not test sensitivity to observed interaction weights, and so
# cannot show that binarisation is inconsequential.
#
# How sensitivity_one_graph() (01_functions.R) works:
#   - the largest connected component (LCC) of the binarised network is taken
#     first: this is the graph classified in the main pipeline;
#   - thresholds are the 0-90% quantiles of the edge overlap values; tied
#     quantiles are dropped, so a network has no row at a level where its
#     overlap values are tied (= not tested at that level);
#   - edges below the threshold are deleted and isolates removed;
#   - "largest_component": the LCC of the pruned graph is re-extracted and
#     classified (the main analysis);
#   - "full_network": the whole pruned graph (from the original LCC, not the
#     original multi-component network) is classified. classify_graph()
#     requires a connected graph, so this mode returns NA once the pruned
#     graph disconnects. It is reported for coverage only.
#
# Networks tested and classifiable, conditional
# label stability, actual % edges removed, fragmentation, by level/class/taxon.
#
# Inputs : data/Networks-Threshold-Test/*.edges, recovered main model
# Outputs: outputs/threshold_sensitivity_{raw,by_graph,wide}.csv   (step 1)
#          outputs/table3_point1_stability_tested.csv              (main Table 3)
#          outputs/tableS1_point1_stability_by_taxon.csv           (Table S1)
#          outputs/tableS_point1_coverage_by_level.csv             (Table S23)
#          outputs/tableS_point1_coverage_by_class.csv             (Table S24)
#          outputs/point1_disconnected_by_level.csv
#          outputs/point1_lcc_extraction.csv
#
#   source(here::here("R", "06_threshold_sensitivity.R"))
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))
set.seed(123)

## Set to TRUE to re-run the (slow) edge-removal test; otherwise the existing
## outputs/threshold_sensitivity_raw.csv is summarised.
RERUN    <- FALSE
RAW_FILE <- here::here("outputs", "threshold_sensitivity_raw.csv")

hr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
LEVELS <- paste0(seq(0, 90, 10), "%")

# =============================================================================
# (1) Edge-removal test
# =============================================================================
if (RERUN || !file.exists(RAW_FILE)) {
  out <- threshold_sensitivity_test(
    edges_dir     = here::here("data", "Networks-Threshold-Test"),
    model         = load_main_model(),
    weight_method = "adjacency",   # neighbourhood overlap (Jaccard)
    n_thresholds  = 10,
    out_prefix    = here::here("outputs", "threshold_sensitivity")
  )
}
raw <- utils::read.csv(RAW_FILE, stringsAsFactors = FALSE)

graphs <- raw %>% dplyr::distinct(graph_name, baseline_class)
cat("Networks:", nrow(graphs), "\n"); print(table(graphs$baseline_class))

## One row per network x level, both modes side by side
w <- raw %>%
  dplyr::select(graph_name, baseline_class, component_mode, threshold_quant,
                n_nodes, n_edges, predicted_class) %>%
  tidyr::pivot_wider(names_from = component_mode,
                     values_from = c(n_nodes, n_edges, predicted_class))

base <- w %>% dplyr::filter(threshold_quant == "0%") %>%
  dplyr::select(graph_name, nodes0 = n_nodes_largest_component, edges0 = n_edges_full_network)

w <- w %>% dplyr::left_join(base, by = "graph_name") %>%
  dplyr::mutate(
    pct_edges_removed = 100 * (1 - n_edges_full_network / edges0),
    # fragmented: the LCC of the pruned graph no longer holds all original nodes
    fragmented        = n_nodes_largest_component < nodes0,
    # disconnected: the whole pruned graph has more than one component
    disconnected      = n_nodes_largest_component < n_nodes_full_network,
    pct_nodes_in_lcc  = 100 * n_nodes_largest_component / nodes0,
    classifiable      = !is.na(predicted_class_largest_component),
    same_class        = classifiable & predicted_class_largest_component == baseline_class,
    full_classifiable = !is.na(predicted_class_full_network),
    taxon = dplyr::case_when(
      grepl("^insecta", graph_name, ignore.case = TRUE) ~ "Ant / Insecta",
      grepl("^chimp",   graph_name, ignore.case = TRUE) ~ "Chimpanzee",
      grepl("^badger",  graph_name, ignore.case = TRUE) ~ "Badger",
      grepl("^giraffe", graph_name, ignore.case = TRUE) ~ "Giraffe",
      grepl("^aves",    graph_name, ignore.case = TRUE) ~ "Aves",
      TRUE ~ "Mammalia (other)"))

# =============================================================================
# (2) "full_network" mode at the 50% label, initially geometric networks
# =============================================================================
hr("(2) 'full_network' mode at the 50% label, initially geometric networks")
geo <- dplyr::filter(graphs, baseline_class == "Spatial")
g50 <- dplyr::filter(w, baseline_class == "Spatial", threshold_quant == "50%")
cat("Initially geometric networks:", nrow(geo), "\n")
cat("  with a row at 50%:          ", nrow(g50), "  (no row:", nrow(geo) - nrow(g50), ")\n")
cat("  full_network unclassifiable:", sum(!g50$full_classifiable), "\n")
cat("  ...of which disconnected:   ", sum(!g50$full_classifiable & g50$disconnected), "\n")

disc_level <- w %>%
  dplyr::mutate(initially_geometric = baseline_class == "Spatial") %>%
  dplyr::group_by(initially_geometric, threshold_quant) %>%
  dplyr::summarise(tested = dplyr::n(), n_disconnected = sum(disconnected),
                   full_classifiable = sum(full_classifiable), .groups = "drop")
print(as.data.frame(disc_level), row.names = FALSE)

# =============================================================================
# (3) Coverage by nominal level, largest-component mode (Table S23)
# =============================================================================
cov_level <- function(ids, label) {
  w %>% dplyr::filter(graph_name %in% ids) %>%
    dplyr::group_by(threshold_quant) %>%
    dplyr::summarise(
      group = label,
      n_networks = length(unique(ids)),
      n_with_row = dplyr::n(),
      n_classifiable = sum(classifiable),
      pct_same_class_of_classifiable = round(100 * sum(same_class) / sum(classifiable), 1),
      median_pct_edges_removed = round(stats::median(pct_edges_removed), 1),
      iqr_pct_edges_removed = paste0(round(stats::quantile(pct_edges_removed, .25), 1), "-",
                                     round(stats::quantile(pct_edges_removed, .75), 1)),
      pct_fragmented = round(100 * mean(fragmented), 1),
      median_pct_nodes_in_lcc = round(stats::median(pct_nodes_in_lcc), 1),
      .groups = "drop") %>%
    dplyr::mutate(threshold_quant = factor(threshold_quant, LEVELS)) %>%
    dplyr::arrange(threshold_quant) %>%
    dplyr::relocate(group)
}
tab_level <- dplyr::bind_rows(cov_level(graphs$graph_name, "All networks"),
                              cov_level(geo$graph_name, "Initially geometric"))
hr("(3) Coverage by nominal removal level (Table S23)")
print(as.data.frame(tab_level), row.names = FALSE)

# =============================================================================
# (4) Coverage by initial class at the 30% and 50% labels (Table S24)
# =============================================================================
by_class <- w %>%
  dplyr::filter(threshold_quant %in% c("30%", "50%")) %>%
  dplyr::group_by(baseline_class, threshold_quant) %>%
  dplyr::summarise(
    n_with_row = dplyr::n(),
    n_classifiable = sum(classifiable),
    pct_same_class = round(100 * sum(same_class) / sum(classifiable), 1),
    median_pct_edges_removed = round(stats::median(pct_edges_removed), 1),
    pct_fragmented = round(100 * mean(fragmented), 1),
    .groups = "drop") %>%
  dplyr::left_join(dplyr::count(graphs, baseline_class, name = "n_networks"), by = "baseline_class") %>%
  dplyr::relocate(n_networks, .after = threshold_quant)
hr("(4) By initial class at the 30% and 50% labels (Table S24)")
print(as.data.frame(by_class), row.names = FALSE)

# =============================================================================
# (5) Label stability among networks tested at 30% and 50%:
#     main Table 3 (by initial class) and Table S1 (by taxon)
# =============================================================================
stab <- function(d, grp) {
  n0 <- d %>% dplyr::filter(threshold_quant == "0%") %>%
    dplyr::group_by(.data[[grp]]) %>%
    dplyr::summarise(N = dplyr::n(), median_nodes = stats::median(nodes0), .groups = "drop")
  lv <- d %>% dplyr::filter(threshold_quant %in% c("30%", "50%")) %>%
    dplyr::group_by(.data[[grp]], threshold_quant) %>%
    dplyr::summarise(tested = dplyr::n(), classifiable = sum(classifiable),
                     pct_same = round(100 * sum(same_class) / sum(classifiable), 1),
                     n_disconnected = sum(disconnected), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = threshold_quant,
                       values_from = c(tested, classifiable, pct_same, n_disconnected))
  dplyr::left_join(n0, lv, by = grp)
}
w_all   <- dplyr::mutate(w, all = "Overall")
table3  <- dplyr::bind_rows(stab(w, "baseline_class"),
                            dplyr::rename(stab(w_all, "all"), baseline_class = all))
tableS1 <- dplyr::bind_rows(stab(w, "taxon"),
                            dplyr::rename(stab(w_all, "all"), taxon = all))
hr("(5) Main Table 3 (by initial class) and Table S1 (by taxon)")
print(as.data.frame(table3), row.names = FALSE)
cat("\n"); print(as.data.frame(tableS1), row.names = FALSE)

# =============================================================================
# (6) Original networks: components and share of nodes in the LCC
# =============================================================================
hr("(6) Original networks: components and share of nodes in the largest component")
files <- list.files(here::here("data", "Networks-Threshold-Test"), pattern = "\\.edges$", full.names = TRUE)
files <- files[tools::file_path_sans_ext(basename(files)) %in% graphs$graph_name]
lcc_tab <- dplyr::bind_rows(lapply(files, function(f) {
  g <- read_edges_graph(f); cs <- igraph::components(g)$csize
  tibble::tibble(graph_name = tools::file_path_sans_ext(basename(f)),
                 n_nodes = igraph::vcount(g), n_components = length(cs),
                 pct_nodes_in_lcc = 100 * max(cs) / sum(cs))
}))
cat("Networks checked:", nrow(lcc_tab), "\n")
cat("Connected (one component):", sum(lcc_tab$n_components == 1), "\n")
cat("Median % of nodes in the largest component:", round(stats::median(lcc_tab$pct_nodes_in_lcc), 1),
    " | 10th percentile:", round(stats::quantile(lcc_tab$pct_nodes_in_lcc, .1), 1), "\n")

# =============================================================================
# Write outputs
# =============================================================================
utils::write.csv(table3,     here::here("outputs", "table3_point1_stability_tested.csv"),  row.names = FALSE)
utils::write.csv(tableS1,    here::here("outputs", "tableS1_point1_stability_by_taxon.csv"), row.names = FALSE)
utils::write.csv(tab_level,  here::here("outputs", "tableS_point1_coverage_by_level.csv"), row.names = FALSE)
utils::write.csv(by_class,   here::here("outputs", "tableS_point1_coverage_by_class.csv"), row.names = FALSE)
utils::write.csv(disc_level, here::here("outputs", "point1_disconnected_by_level.csv"),    row.names = FALSE)
utils::write.csv(lcc_tab,    here::here("outputs", "point1_lcc_extraction.csv"),           row.names = FALSE)

message("06_threshold_sensitivity.R complete.")