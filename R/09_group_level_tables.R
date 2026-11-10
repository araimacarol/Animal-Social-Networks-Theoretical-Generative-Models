# =============================================================================
# 09_group_level_tables.R
# -----------------------------------------------------------------------------
# Writes the contents of the group-level tables: main Table 4 and Tables S12
# (badgers), S13 (chimpanzees) and S14 (ants, aggregate and daily
# networks), all with the RECOVERED first-stage classifier.
#
# For each network:
#   - number of nodes in the largest connected component (the graph analysed);
#   - modularity of the whole network (walktrap partition, as in script 05);
#   - predicted class of the whole network: taken from script 02's output
#     (outputs/Predicted-Class-of-Animal-Social-network.csv), the same
#     classification used everywhere else in the paper; if a network is not in
#     that file, the class computed here from the edge list is used and flagged;
#   - every group of at least MIN_GROUP_SIZE nodes and its predicted class.
#     Groups: field-defined (prescribed) groups for badgers, walktrap for
#     chimpanzees, giraffes and ants -- the same rules as the manuscript.
#
# A single walktrap community spanning the whole network is not a subgroup; it
# is written with is_whole_network = TRUE and shown as "None" in the tables.
#
# Inputs : data/Empirical-Network-Edges/*.edges, data/Badger_Groups.csv,
#          outputs/Predicted-Class-of-Animal-Social-network.csv (script 02)
# Outputs: outputs/group_tables_networks.csv   (one row per network)
#          outputs/group_tables_groups.csv     (one row per classified group)
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

model <- load_main_model()
LABEL <- c(ER = "Erdos-Renyi", sbm = "SBM", SF = "Scale-free",
           Spatial = "Geometric graph", SW = "Small-world")

## ---- Networks shown in the tables -------------------------------------------
BADGERS <- c(paste0("Badger_Spring", 1:3), paste0("Badger_Summer", 1:4),
             paste0("Badger_Autumn", 1:4), paste0("Badger_Winter", 2:4))
CHIMPS  <- c(paste0("Chimp_", c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Dec"),
                    "_Prob5_MatrixNoNAs"),                                       # 5 m
             paste0("Chimp_", c("Jan", "Feb", "Mar", "Apr", "May"), "SRI_G_Matrix"),
             paste0("Chimp__", c("Jun", "Jul", "Aug", "Dec"), "SRI_G_Matrix"))   # 50 m
CHIMPS  <- sub("^Chimp__Jul", "Chimp_Jul", CHIMPS)   # file is Chimp_JulSRI_G_Matrix
GIRAFFE <- "Giraffe_edgelist"
ANTS    <- unlist(lapply(1:6, function(k) c(
  sprintf("insecta-ant-colony%d", k),
  sprintf("insecta-ant-colony%d-day%s", k, c("05", "20", "35")))))
ALL <- c(BADGERS, CHIMPS, GIRAFFE, ANTS)

paths   <- vapply(ALL, resolve_edge_file, character(1), dir = dir_edges)
missing <- ALL[is.na(paths)]
if (length(missing) > 0)
  message("No edge file for (skipped): ", paste(missing, collapse = ", "))
paths <- paths[!is.na(paths)]

clean_name <- function(x) sub("\\.csv$", "", trimws(x))

## ---- Groups (walktrap is deterministic, so no seed is needed) ---------------
cd <- summarize_group_sizes_and_predictions(
  edge_files = unname(paths), model = model,
  methods = c("walktrap", "prescribed"), min_group_size = MIN_GROUP_SIZE)

groups <- cd$predictions %>%
  dplyr::mutate(graph_name = clean_name(graph_name)) %>%
  dplyr::filter(ifelse(grepl("^Badger", graph_name), method == "prescribed",
                       method == "walktrap"))

## ---- Whole-network size, modularity and class -------------------------------
stage1_file <- here::here("outputs", "Predicted-Class-of-Animal-Social-network.csv")
stage1 <- if (file.exists(stage1_file)) utils::read.csv(stage1_file) else
  data.frame(target = character(), Predicted_classes = character())
stage1$key <- gsub("[^a-z0-9]", "", tolower(clean_name(stage1$target)))

networks <- dplyr::bind_rows(lapply(names(paths), function(nm) {
  df <- utils::read.table(paths[[nm]], sep = " ", fill = TRUE, header = FALSE)
  g  <- igraph::graph_from_data_frame(as.matrix(df), directed = FALSE)
  g  <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g  <- largest_component(igraph::delete_vertices(g, which(igraph::degree(g) == 0)))
  key <- gsub("[^a-z0-9]", "", tolower(nm))
  s1  <- stage1$Predicted_classes[stage1$key == key]
  own <- classify_graph(g, model, tag = nm)
  tibble::tibble(
    graph_name = nm, n_nodes = igraph::vcount(g),
    modularity = round(igraph::modularity(igraph::cluster_walktrap(g)), 4),
    class_stage1 = if (length(s1) == 1) as.character(s1) else NA_character_,
    class_from_edges = if (is.null(own)) NA_character_ else own)
})) %>%
  dplyr::mutate(
    predicted_class = dplyr::coalesce(class_stage1, class_from_edges),
    note = dplyr::case_when(
      is.na(class_stage1) ~ "not in script-02 output; class computed from edge list",
      class_stage1 != class_from_edges ~ "edge-list features give a different class",
      TRUE ~ ""),
    predicted_label = LABEL[predicted_class])

groups <- groups %>%
  dplyr::left_join(dplyr::select(networks, graph_name, n_nodes_network = n_nodes),
                   by = "graph_name") %>%
  dplyr::mutate(is_whole_network = n_nodes == n_nodes_network,
                predicted_label = LABEL[predicted_class]) %>%
  dplyr::arrange(graph_name, dplyr::desc(n_nodes))

## ---- Checks and output ---------------------------------------------------------
cat("\nNetworks:", nrow(networks), "| classified groups:", sum(!groups$is_whole_network), "\n")
if (any(networks$note != "")) {
  cat("\nNetworks needing attention:\n")
  print(as.data.frame(dplyr::filter(networks, note != "")), row.names = FALSE)
}
print(as.data.frame(networks), row.names = FALSE)
print(as.data.frame(groups), row.names = FALSE)

utils::write.csv(networks, here::here("outputs", "group_tables_networks.csv"), row.names = FALSE)
utils::write.csv(groups,   here::here("outputs", "group_tables_groups.csv"),   row.names = FALSE)

message("09_group_level_tables.R complete.")
