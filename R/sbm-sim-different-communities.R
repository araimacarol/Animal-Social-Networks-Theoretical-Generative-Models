# =============================================================================
# sbm-sim-different-communities.R  --  SBM simulations (Supp. Tables S16/S17)
# -----------------------------------------------------------------------------
# Simulates stochastic block models (SBMs) with k = 2-6 communities and
# n = 100-1000 nodes, under lower and higher between-community connection
# probability, and classifies each network with the first-stage classifier.
#
# Design (unchanged from the run reported in the paper):
#   - one network per k and n (n = 100, 150, ..., 1000)
#   - block sizes floor(n/k), remainder added to block 1
#   - within-community probability p ~ U(0.5, 0.8), drawn per community
#   - between-community probability q drawn per pair of communities:
#       higher: q ~ U(0.1, 0.4)    lower: q ~ U(0.01, 0.04)
#   - only connected networks are classified
#
# Outputs (one row per network; features and predicted class from the same
# network):
#   outputs/all_sbm_preds_higher_between_connectn_probs.csv  -> Table S17
#   outputs/all_sbm_preds_lower_between_connectn_probs.csv   -> Table S16
# The tables are written from these files by 08_sbm_simulation_tables.R.
#
# Note: the run reported in the paper did not set an initial seed, so a rerun
# produces new networks from the same design rather than the identical ones.
# The deposited CSV files above are those used in the paper; to protect them,
# this script writes a rerun to the same names with the suffix "_rerun".
#
# Not part of run_all.R. Run on its own (slow for large n):
#   source(here::here("R", "sbm-sim-different-communities.R"))
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

set.seed(1)
main_model <- load_main_model()

NODE_SIZES  <- seq(100, 1000, by = 50)
K_VALUES    <- 2:6
N_INSTANCES <- 1

## ---- Probability matrix: p within communities, q between each pair ---------
sbm_prob_matrix <- function(k, p_between_range, p_within_range = c(0.5, 0.8)) {
  mat <- matrix(0, nrow = k, ncol = k)
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      mat[i, j] <- mat[j, i] <- stats::runif(1, p_between_range[1], p_between_range[2])
    }
  }
  diag(mat) <- stats::runif(k, p_within_range[1], p_within_range[2])
  mat
}

## ---- One SBM with (near-)equal block sizes ----------------------------------
make_sbm <- function(pm, n, k) {
  block_sizes <- rep(floor(n / k), k)
  block_sizes[1] <- block_sizes[1] + (n - sum(block_sizes))
  g <- igraph::sample_sbm(n, pref.matrix = pm, block.sizes = block_sizes)
  igraph::V(g)$community <- rep(seq_len(k), times = block_sizes)
  g
}

## ---- Simulate, compute features, classify -----------------------------------
generate_sbm_networks <- function(node_sizes, k_values, n_instances, p_between_range) {
  out <- list()
  for (n in node_sizes) {
    for (k in k_values) {
      for (i in seq_len(n_instances)) {
        pm <- sbm_prob_matrix(k, p_between_range)
        stopifnot(isSymmetric(pm))
        g  <- make_sbm(pm, n, k)
        message(sprintf("k = %d, n = %d, instance %d", k, n, i))
        
        feats <- calc_graph_features(list(g))
        if (feats$connected != 1) next          # classify connected networks only
        feats$graph_name <- "sbm"
        
        sim_data <- feats %>%
          dplyr::select(dplyr::all_of(FEATURE_COLS)) %>%
          dplyr::mutate_if(is.character, factor) %>%
          janitor::clean_names() %>%
          tibble::as_tibble()
        
        out[[length(out) + 1]] <- tibble::tibble(
          k_communities = k,
          .pred_class   = predict_main_model(main_model, sim_data)$predicted_class) %>%
          dplyr::bind_cols(sim_data) %>%
          dplyr::mutate(n_nodes = n, instance = i,
                        p_within = paste(diag(pm), collapse = ","))
      }
    }
  }
  dplyr::bind_rows(out)
}

## ---- Higher between-community probability (Table S17) -----------------------
preds_higher <- generate_sbm_networks(NODE_SIZES, K_VALUES, N_INSTANCES,
                                      p_between_range = c(0.1, 0.4))
utils::write.csv(preds_higher,
                 here::here("outputs", "all_sbm_preds_higher_between_connectn_probs_rerun.csv"),
                 row.names = FALSE)
cat("Higher q: mean modularity", round(mean(preds_higher$modularity), 3), "\n")
print(table(preds_higher$.pred_class))

## ---- Lower between-community probability (Table S16) ------------------------
preds_lower <- generate_sbm_networks(NODE_SIZES, K_VALUES, N_INSTANCES,
                                     p_between_range = c(0.01, 0.04))
utils::write.csv(preds_lower,
                 here::here("outputs", "all_sbm_preds_lower_between_connectn_probs_rerun.csv"),
                 row.names = FALSE)
cat("Lower q: mean modularity", round(mean(preds_lower$modularity), 3), "\n")
print(table(preds_lower$.pred_class))

message("sbm-sim-different-communities.R complete.")