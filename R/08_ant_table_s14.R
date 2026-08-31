# =============================================================================
# 08_ant_table_s14.R
# -----------------------------------------------------------------------------
# Rebuilds Table S14: the ant (C. fellah) networks and their groups -- main
# network modularity and predicted class, plus the predicted class of every
# group of at least eight nodes.
#
# WHY THIS SCRIPT EXISTS
#   The published Table S14 was produced before the first-stage model was
#   recovered to a version-portable JSON booster, and before the igraph 2.0
#   return-type fix in classify_subgraphs(). It contained arithmetic
#   impossibilities (a 97-node group inside a 73-node network) and rows that
#   disagreed with Table S6. This script regenerates the whole table from the
#   same recovered model used by scripts 02, 05 and 06.
#
#   Walktrap only -- Table S14 reports one partitioning method per network. The
#   cross-algorithm comparison lives in script 05 / Table S6.
#
# Inputs : data/Empirical-Network-Edges/insecta-ant-colony*.edges,
#          recovered first-stage model
# Outputs: outputs/table_s14_ants.csv
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

set.seed(123)
model <- load_main_model()

## ---- The 24 ant networks in Table S14 ---------------------------------------
## Six aggregate colony networks plus days 05, 20 and 35 for each colony.
ANT_NETWORKS <- c(
  paste0("insecta-ant-colony", 1:6),
  paste0("insecta-ant-colony", rep(1:6, each = 3), "-day",
         rep(c("05", "20", "35"), times = 6))
)

ant_files <- resolve_edge_files(ANT_NETWORKS, dir = dir_edges)
stopifnot(length(ant_files) == length(ANT_NETWORKS))
message("Rebuilding Table S14 from ", length(ant_files), " ant networks.")

## ---- One network: modularity, main class, group classes ---------------------
s14_one_network <- function(file, model, min_group_size = MIN_GROUP_SIZE) {
  gname <- tools::file_path_sans_ext(basename(file))
  message("Processing: ", gname)

  g <- tryCatch(read_edges_graph(file), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  g <- largest_component(g)
  if (igraph::vcount(g) < min_group_size) return(NULL)

  cw         <- igraph::cluster_walktrap(g)
  membership <- igraph::membership(cw)
  mod_val    <- igraph::modularity(g, membership)
  main_class <- classify_graph(g, model, tag = gname)
  if (is.null(main_class)) main_class <- NA_character_

  grp_rows <- list()
  for (grp in unique(membership)) {
    nodes <- which(membership == grp)
    if (length(nodes) < min_group_size) next
    sub <- igraph::induced_subgraph(g, nodes)
    if (!igraph::is_connected(sub) || igraph::ecount(sub) == 0) next
    p <- classify_graph(sub, model, tag = paste0(gname, "_g", grp))
    grp_rows[[length(grp_rows) + 1]] <- tibble::tibble(
      group_id         = as.character(grp),
      group_nodes      = as.integer(igraph::vcount(sub)),
      group_prediction = if (is.null(p)) NA_character_ else p)
  }

  ## A single community spanning the whole network is reported as one group,
  ## not as "no groups" -- it clears the size threshold.
  grp <- if (length(grp_rows) == 0)
    tibble::tibble(group_id = NA_character_, group_nodes = NA_integer_,
                   group_prediction = paste0("None (>= ", min_group_size, ")"))
  else dplyr::bind_rows(grp_rows)

  tibble::tibble(network         = gname,
                 n_nodes         = as.integer(igraph::vcount(g)),
                 modularity      = round(mod_val, 4),
                 main_prediction = main_class) %>%
    dplyr::slice(rep(1, nrow(grp))) %>%
    dplyr::bind_cols(grp)
}

## ---- Run --------------------------------------------------------------------
s14 <- dplyr::bind_rows(lapply(ant_files, s14_one_network, model = model))

relabel <- c("ER" = "Erdos-Renyi", "sbm" = "SBM", "SF" = "Scale Free",
             "Spatial" = "Geometric Graph", "SW" = "Small World")
s14 <- s14 %>%
  dplyr::mutate(main_prediction  = dplyr::recode(main_prediction,  !!!relabel),
                group_prediction = dplyr::recode(group_prediction, !!!relabel))

## ---- Internal consistency checks -------------------------------------------
## The published table had groups larger than their own network. Fail loudly
## rather than shipping that again.
bad_size <- s14 %>%
  dplyr::filter(!is.na(group_nodes), group_nodes > n_nodes)
if (nrow(bad_size) > 0) {
  print(as.data.frame(bad_size), row.names = FALSE)
  stop("Group larger than its network -- check the partitioning.")
}

bad_sum <- s14 %>%
  dplyr::filter(!is.na(group_nodes)) %>%
  dplyr::group_by(network, n_nodes) %>%
  dplyr::summarise(sum_groups = sum(group_nodes), .groups = "drop") %>%
  dplyr::filter(sum_groups > n_nodes)
if (nrow(bad_sum) > 0) {
  print(as.data.frame(bad_sum), row.names = FALSE)
  stop("Group sizes sum to more than the network -- check the partitioning.")
}

message("Consistency checks passed.")

print(as.data.frame(s14), row.names = FALSE)
utils::write.csv(s14, here::here("outputs", "table_s14_ants.csv"),
                 row.names = FALSE, na = "")
message("Wrote: ", here::here("outputs", "table_s14_ants.csv"))

## ---- Summary for the Results text -------------------------------------------
daily <- s14 %>% dplyr::filter(grepl("-day", network))
aggreg <- s14 %>% dplyr::filter(!grepl("-day", network))

cat("\n### AGGREGATE COLONY NETWORKS ###\n")
print(aggreg %>% dplyr::distinct(network, n_nodes, modularity, main_prediction) %>%
        as.data.frame(), row.names = FALSE)

cat("\n### DAILY NETWORKS: main-network class ###\n")
print(daily %>% dplyr::distinct(network, modularity, main_prediction) %>%
        dplyr::count(main_prediction, name = "n_networks") %>% as.data.frame(),
      row.names = FALSE)

cat("\n### DAILY NETWORKS: group-level class ###\n")
print(daily %>% dplyr::filter(!is.na(group_nodes)) %>%
        dplyr::count(group_prediction, name = "n_groups") %>% as.data.frame(),
      row.names = FALSE)

message("08_ant_table_s14.R complete.")
