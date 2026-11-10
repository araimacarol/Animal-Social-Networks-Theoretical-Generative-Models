# =============================================================================
# PCA_Applicability.R
# -----------------------------------------------------------------------------
# Applicability of the first-stage classifier to smaller
# empirical networks. No retraining. Answers:
#
#   (a)  Did the training networks, after preprocessing, cover networks below
#        50 nodes?
#   (b)  If not, validate across the empirical size and density ranges
#        (simulated networks, n = 11-45, original generators).
#   (c)  Two PCs (60.6%) cannot establish coverage -> check all 18 features.
#   (d)  Define the 95.9% coverage criterion (convex hull, Figure 2 PCA).
#   (e)  Classification uncertainty for the empirical networks.
#
# Run from the project (it sources 00_setup.R and 01_functions.R itself):
#   source(here::here("R", "PCA_Applicability.R"))
#
# Inputs (project-relative):
#   Parts (a), (c), (d) need the TRAINING data (a model file holds no training
#   data, so the recovered booster cannot replace them). Two routes:
#   - small files, written once by point2_export_derived.R (preferred):
#       models/point2_training_size_by_class.csv    (a)
#       models/point2_training_feature_ranges.csv   (c)
#       models/point2_fig2_pca.rds                  (d)  Figure 2, Fig. S1
#   - otherwise the original large files:
#       models/df_train_new.rds         309,499 pre-SMOTE training networks (a, c)
#       models/final.fitted.model.rds   fitted workflow; only its stored training
#                                       predictors (mold) are read, its booster is
#                                       never used for prediction             (d)
#   models/final.emp.data.rds            636 empirical networks                (a, c, d)
#   models/booster_portable_mainmodel.json + preproc_trained_mainmodel.rds
#                                        recovered classifier: ALL predictions (b, e)
#   data/GraphFeatOnAllAnimNets.csv      636 networks after order > 10         (e)
#
# Outputs (outputs/ and figures/):
#   tableS_point2a_training_size_range.csv
#   tableS_point2b_confusion.csv, tableS_point2b_accuracy_by_size.csv,
#   tableS_point2b_accuracy_by_density.csv
#   tableS_point2c_feature_range.csv, tableS_point2c_all_feature_summary.csv
#   point2d_pca_convex_hull_coverage.csv
#   tableS_point2e_uncertainty_summary.csv
#   dataS_point2e_empirical_class_probabilities.csv
#   point2_validation_small_networks.csv, point2_validation_draws.csv
#   figures/pca_emp_vs_sim_new.png (Figure 2), figures/pca_feature_loadings.png (Fig. S1)
#
# Manuscript: Table S18 = (a), S19-S20 = (b), S21 = (c), Figure 2 / Fig. S1 /
#   95.9% = (d), S22 and S25 = (e).
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

DATA_DIR <- here::here("models")
OUT_DIR  <- here::here("outputs")
FIG_DIR  <- here::here("figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

hr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

emp_fig <- readRDS(file.path(DATA_DIR, "final.emp.data.rds"))
DERIVED <- file.path(DATA_DIR, c("point2_training_size_by_class.csv",
                                 "point2_training_feature_ranges.csv",
                                 "point2_fig2_pca.rds"))
USE_DERIVED <- all(file.exists(DERIVED))
if (!USE_DERIVED) {
  need <- file.path(DATA_DIR, c("df_train_new.rds", "final.fitted.model.rds"))
  if (!all(file.exists(need)))
    stop("Run point2_export_derived.R once, or put these in models/:\n  ",
         paste(basename(need[!file.exists(need)]), collapse = "\n  "))
  train_df <- readRDS(file.path(DATA_DIR, "df_train_new.rds"))
}
cat("Empirical networks:", nrow(emp_fig), "| training data from",
    if (USE_DERIVED) "derived files" else "df_train_new.rds / final.fitted.model.rds", "\n")


# =============================================================================
# (a) TRAINING SIZE RANGE
# =============================================================================
hr("(a) Did the training networks cover networks below 50 nodes?")

if (USE_DERIVED) {
  size_tab <- utils::read.csv(DERIVED[1])
} else {
  cls_col  <- intersect(c("graph_name", "target", "GraphName"), names(train_df))[1]
  size_tab <- train_df %>%
    dplyr::group_by(class = as.character(.data[[cls_col]])) %>%
    dplyr::summarise(n_networks = dplyr::n(), min_order = min(order), max_order = max(order),
                     n_order_lt_50 = sum(order < 50), .groups = "drop") %>%
    dplyr::bind_rows(dplyr::summarise(train_df, class = "All", n_networks = dplyr::n(),
                                      min_order = min(order), max_order = max(order),
                                      n_order_lt_50 = sum(order < 50)))
}
print(size_tab)
n_emp_lt50 <- sum(emp_fig$order < 50)
cat("Empirical networks below 50 nodes:", n_emp_lt50, "of", nrow(emp_fig), "\n")
utils::write.csv(size_tab, file.path(OUT_DIR, "tableS_point2a_training_size_range.csv"), row.names = FALSE)


# =============================================================================
# (c) FULL 18-FEATURE RANGE CHECK
# =============================================================================
hr("(c) Coverage in the full feature space: empirical vs training range, all features")

if (USE_DERIVED) {
  rng <- utils::read.csv(DERIVED[2])
  feat_names <- rng$feature
  lo <- stats::setNames(rng$train_min, feat_names); hi <- stats::setNames(rng$train_max, feat_names)
} else {
  feat_names <- intersect(setdiff(names(emp_fig), "graph_name"), names(train_df))
  lo     <- sapply(train_df[feat_names], min, na.rm = TRUE)
  hi     <- sapply(train_df[feat_names], max, na.rm = TRUE)
}
inside <- sapply(feat_names, function(f) emp_fig[[f]] >= lo[f] & emp_fig[[f]] <= hi[f])
all_in <- rowSums(inside) == length(feat_names)

range_tab <- tibble::tibble(feature = feat_names, train_min = lo, train_max = hi,
                            n_empirical_outside = colSums(!inside)) %>%
  dplyr::arrange(dplyr::desc(n_empirical_outside))
band_fig <- cut(emp_fig$order, c(10, 49, 1000, Inf), labels = c("<50", "50-1000", ">1000"))
all_tab  <- tibble::tibble(group = c("All", levels(band_fig)),
                           n = c(nrow(emp_fig), as.vector(table(band_fig))),
                           n_inside_all_features = c(sum(all_in), tapply(all_in, band_fig, sum)))

cat("Features checked:", length(feat_names), "\n\nOutside the training range, per feature:\n")
print(range_tab, n = Inf)
cat("\nInside the training range on all features:\n"); print(all_tab)
utils::write.csv(range_tab, file.path(OUT_DIR, "tableS_point2c_feature_range.csv"), row.names = FALSE)
utils::write.csv(all_tab,   file.path(OUT_DIR, "tableS_point2c_all_feature_summary.csv"), row.names = FALSE)

if (!USE_DERIVED) rm(train_df)


# =============================================================================
# (d) FIGURE 2 PCA, CONVEX-HULL COVERAGE, FIGURES
# -----------------------------------------------------------------------------
# Same data and PCA as the published Figure 2 (workflow mold + empirical,
# normalise, 2 components). Published: PC1 = 38.3%, PC2 = 22.3%, hull 95.9%.
# =============================================================================
hr("(d) Definition of the 95.9%: convex hull in the Figure 2 PCA")

if (USE_DERIVED) {
  pca_out <- readRDS(DERIVED[3])
  pc1_var <- round(pca_out$variance[["PC1"]], 1); pc2_var <- round(pca_out$variance[["PC2"]], 1)
  sim_pts <- dplyr::mutate(tibble::as_tibble(pca_out$sim), Data_Source = "Classifier Training")
  emp_pts <- dplyr::mutate(dplyr::select(pca_out$emp, PC1, PC2), Data_Source = "ASNR Benchmark")
  pca_loadings_raw <- pca_out$loadings
} else {
  fit  <- readRDS(file.path(DATA_DIR, "final.fitted.model.rds"))
  wf   <- tryCatch(tune::extract_workflow(fit), error = function(e) fit$.workflow[[1]])
  mold <- hardhat::extract_mold(wf)$predictors %>% tibble::as_tibble()
  final_features <- colnames(mold)
  rm(fit, wf)
  combined_data <- dplyr::bind_rows(
    dplyr::mutate(mold, Data_Source = "Classifier Training"),
    emp_fig %>% dplyr::select(dplyr::all_of(final_features)) %>%
      dplyr::mutate(Data_Source = "ASNR Benchmark"))
  rm(mold)
  pca_prep <- recipes::recipe(~ ., data = combined_data) %>%
    recipes::update_role(Data_Source, new_role = "id") %>%
    recipes::step_normalize(recipes::all_predictors()) %>%
    recipes::step_pca(recipes::all_predictors(), num_comp = 2) %>%
    recipes::prep()
  pca_baked <- recipes::bake(pca_prep, new_data = NULL)
  pca_variance <- recipes::tidy(pca_prep, number = 2, type = "variance") %>%
    dplyr::filter(terms == "percent variance")
  pc1_var <- round(pca_variance$value[pca_variance$component == 1], 1)
  pc2_var <- round(pca_variance$value[pca_variance$component == 2], 1)
  sim_pts <- dplyr::filter(pca_baked, Data_Source == "Classifier Training")
  emp_pts <- dplyr::filter(pca_baked, Data_Source == "ASNR Benchmark")
  pca_loadings_raw <- recipes::tidy(pca_prep, number = 2, type = "coef") %>%
    dplyr::filter(component %in% c("PC1", "PC2")) %>% dplyr::select(terms, value, component)
  rm(combined_data, pca_baked)
}

hull    <- sim_pts[grDevices::chull(sim_pts$PC1, sim_pts$PC2), ]
in_hull <- sp::point.in.polygon(emp_pts$PC1, emp_pts$PC2, hull$PC1, hull$PC2) > 0

cov_tab <- tibble::tibble(
  pc1_pct = pc1_var, pc2_pct = pc2_var,
  n_simulated = nrow(sim_pts), n_empirical = nrow(emp_pts),
  n_inside_convex_hull = sum(in_hull), pct_inside_convex_hull = round(100 * mean(in_hull), 1))
print(cov_tab)
utils::write.csv(cov_tab, file.path(OUT_DIR, "point2d_pca_convex_hull_coverage.csv"), row.names = FALSE)

## ---- Figure 2: empirical vs simulated --------------------------------------
COL_SIM <- "#9ecae1"; COL_DEN <- "#377eb8"; COL_EMP <- "#e41a1c"
LAB <- c("Classifier Training" = "Simulated (Trained)", "ASNR Benchmark" = "Empirical ASNs")

## Simulated points are drawn first and the empirical networks on top. (Earlier
## versions also drew 2-D density contours of the training networks; they lay
## beneath the ~474,000 simulated points and were not visible, so they are omitted.)
p_pca <- ggplot2::ggplot(mapping = ggplot2::aes(x = PC1, y = PC2)) +
  ggplot2::geom_point(data = sim_pts, ggplot2::aes(colour = Data_Source, shape = Data_Source),
                      alpha = 0.5, size = 1.2) +
  ggplot2::geom_point(data = emp_pts, ggplot2::aes(colour = Data_Source, shape = Data_Source),
                      alpha = 0.5, size = 1.6) +
  ggplot2::scale_colour_manual(name = "Data", labels = LAB,
                               values = c("Classifier Training" = COL_SIM, "ASNR Benchmark" = COL_EMP)) +
  ggplot2::scale_shape_manual(name = "Data", labels = LAB,
                              values = c("Classifier Training" = 16, "ASNR Benchmark" = 17)) +
  ggplot2::guides(colour = ggplot2::guide_legend(title = "Data", override.aes = list(alpha = 1, size = 2.5)),
                  shape = ggplot2::guide_legend(title = "Data")) +
  ggplot2::coord_cartesian(xlim = c(-15, 6), ylim = c(-10, 20)) +
  ggplot2::labs(x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(family = "sans", size = 11),
                 axis.title = ggplot2::element_text(size = 11),
                 axis.text = ggplot2::element_text(size = 10),
                 legend.title = ggplot2::element_text(size = 11),
                 legend.text = ggplot2::element_text(size = 10),
                 legend.background = ggplot2::element_rect(fill = "white", colour = NA),
                 legend.box = "vertical", legend.position = "right")

while (!is.null(grDevices::dev.list())) grDevices::dev.off()
grDevices::png(file.path(FIG_DIR, "pca_emp_vs_sim_new.png"), width = 850, height = 600, res = 120)
print(p_pca); grDevices::dev.off()

## ---- Fig. S1: feature loadings ---------------------------------------------
FEATURE_LABELS <- c(
  "Mean path length" = "mean_path_length", "Mean eccentricity" = "mean_eccentr",
  "Diameter" = "diameter", "Modularity" = "modularity",
  "Eigenvector centrality" = "eigen_centr", "Degree assortativity" = "deg_assort_coef",
  "Betweenness centrality" = "betw_centr", "Order" = "order",
  "Transitivity" = "transitivity", "Closeness centrality" = "closeness_centr",
  "Graph energy" = "graph_energy", "Degree centrality" = "deg_centr",
  "Edges" = "edges", "Mean degree" = "mean_degree", "Min cut" = "min_cut",
  "Fiedler value" = "fiedler_value", "Normalised Fiedler value" = "normalized_fiedler_value",
  "Spectral radius" = "spectral_radius")

pca_loadings <- pca_loadings_raw %>%
  dplyr::mutate(terms = forcats::fct_recode(factor(terms), !!!FEATURE_LABELS),
                terms = tidytext::reorder_within(terms, value, component))

p_load <- ggplot2::ggplot(pca_loadings, ggplot2::aes(x = value, y = terms, fill = value > 0)) +
  ggplot2::geom_col(show.legend = FALSE, width = 0.7) +
  tidytext::scale_y_reordered() +
  ggplot2::facet_wrap(~ component, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(`TRUE` = COL_DEN, `FALSE` = COL_EMP)) +
  ggplot2::labs(x = "Loading coefficient", y = "Network feature") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(family = "sans", size = 11),
                 strip.text = ggplot2::element_text(size = 12, face = "bold"),
                 axis.title = ggplot2::element_text(size = 11),
                 axis.text = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 9.5),
                 panel.spacing = ggplot2::unit(2, "lines"))

grDevices::png(file.path(FIG_DIR, "pca_feature_loadings.png"), width = 1050, height = 650, res = 120)
print(p_load); grDevices::dev.off()

rm(sim_pts, hull)


# =============================================================================
# SHARED: recovered classifier + probability helpers  (used by b and e)
# =============================================================================
model <- load_main_model()

prob_main <- function(model, feature_df) {
  feature_df <- dplyr::mutate(feature_df, dplyr::across(dplyr::where(is.integer), as.double))
  baked <- recipes::bake(model$preproc, new_data = feature_df)
  keep  <- setdiff(names(baked), intersect(c("graph_name", "target", "class", ".pred_class"), names(baked)))
  X <- as.matrix(baked[, keep, drop = FALSE]); storage.mode(X) <- "double"
  P <- stats::predict(model$booster, X)
  if (!is.matrix(P)) P <- matrix(P, nrow = nrow(X), byrow = TRUE)
  colnames(P) <- model$class_order
  P
}

summarise_probs <- function(P) {
  s <- t(apply(P, 1, sort, decreasing = TRUE))
  tibble::tibble(predicted = colnames(P)[max.col(P, ties.method = "first")],
                 p_max = s[, 1], margin = s[, 1] - s[, 2])
}


# =============================================================================
# (b) VALIDATION ON SIMULATED SMALL NETWORKS (n = 11-45)
# -----------------------------------------------------------------------------
# Generators from the classifier-1 data-generation script. Where that script and
# Table S1 disagree, the values matching the stored training data are used:
#   ER       sample_gnp, p = 0.1, ..., 0.9
#   Spatial  fastSpatialNetwork(makeConnected = TRUE), r = 0.1, ..., 0.9
#   SF       sample_pa(psumtree), power in {1, 2}, 1 <= m <= 35
#   SW       sample_smallworld(dim = 1), nei 1-35, p in {0.1, 0.2}
#   SBM      two blocks (0.4n, 0.6n); p_within ~ U(0.5, 0.9), p_between from
#            U(0.1, 0.4) (commented-out values in sbmfunc(), which match the
#            stored training SBM density 0.33-0.65)
# Disconnected draws are redrawn, mirroring filter(connected == 1) in training.
# Generator check (2026-10-10, n = 50 / 100, inside the training domain):
#   SF, SW, Spatial 1.00 / 1.00; ER 0.94 / 0.995; SBM 0.99 / 0.80.
# =============================================================================
hr("(b) Validation on simulated networks below the training floor, n = 11-45")

SIZES     <- c(11, 15, 20, 25, 30, 35, 40, 45)
N_PER     <- 200
MAX_TRIES <- 100
ER_P <- seq(0.1, 0.9, 0.1); SP_R <- seq(0.1, 0.9, 0.1)
SF_POWER <- c(1, 2); SF_M <- 1:35
SW_P <- c(0.1, 0.2);  SW_NEI <- 1:35

## fastSpatialNetwork, copied from the generation script
## (only change: is.connected() -> igraph::is_connected(), same behaviour)
cellCoordsToIndex <- function(i, j, size) (i - 1) * size + j

fastSpatialNetwork <- function(n = 100, r = 0.1, makeConnected = FALSE, keepCellsSeparate = FALSE) {
  v <- data.frame(matrix(nrow = n)); v$x <- runif(n); v$y <- runif(n)
  v[, 1] <- 1:n; colnames(v) <- c("id", "x", "y")
  cellSize <- r / sqrt(2); R2 <- r * r; gridSize <- ceiling(sqrt(2) / r)
  Cell <- vector(mode = "list", gridSize * gridSize)
  for (i in 1:n) {
    idx <- (floor(v$y[i] / cellSize)) * gridSize + floor(v$x[i] / cellSize) + 1
    Cell[[idx]] <- cbind(Cell[[idx]], i)
  }
  E <- matrix(nrow = 0, ncol = 2, byrow = TRUE)
  for (cell in Cell) if (length(cell) > 1) E <- rbind(E, t(combn(cell, 2)))
  proximalPairs <- data.frame(matrix(5, ncol = 3, nrow = 1)); colnames(proximalPairs) <- c("r", "i", "j")
  try_pair <- function(A, B) {
    for (a in A) for (b in B) {
      d <- (v[a, 2] - v[b, 2])^2 + (v[a, 3] - v[b, 3])^2
      if (d < R2) E <<- as.matrix(rbind(E, c(a, b))) else proximalPairs <<- rbind(proximalPairs, c(d, a, b))
    }
  }
  if (keepCellsSeparate == FALSE) {
    for (i in 2:gridSize) {
      try_pair(Cell[[cellCoordsToIndex(1, i, gridSize)]], Cell[[cellCoordsToIndex(1, i - 1, gridSize)]])
      try_pair(Cell[[cellCoordsToIndex(i, 1, gridSize)]], Cell[[cellCoordsToIndex(i - 1, 1, gridSize)]])
      for (j in 2:gridSize) {
        thisIdx <- cellCoordsToIndex(i, j, gridSize); aboveIdx <- cellCoordsToIndex(i, j - 1, gridSize)
        leftIdx <- cellCoordsToIndex(i - 1, j, gridSize); aboveLeftIdx <- cellCoordsToIndex(i - 1, j - 1, gridSize)
        try_pair(Cell[[thisIdx]], Cell[[aboveIdx]])
        try_pair(Cell[[thisIdx]], Cell[[leftIdx]])
        try_pair(Cell[[thisIdx]], Cell[[aboveLeftIdx]])
        try_pair(Cell[[leftIdx]], Cell[[aboveIdx]])
      }
    }
  }
  G <- igraph::make_empty_graph(n, directed = FALSE)
  G <- igraph::add_edges(G, t(E), directed = FALSE)
  if (!makeConnected || igraph::is_connected(G)) return(G)
  proximalPairs <- proximalPairs[-1, ]
  indx <- order(proximalPairs$r)
  comps <- igraph::components(G, mode = "weak")
  for (nu in indx) {
    i <- proximalPairs$i[nu]; j <- proximalPairs$j[nu]
    if (comps$membership[i] != comps$membership[j]) {
      G <- G + igraph::edge(i, j); comps <- igraph::components(G, mode = "weak")
    }
    if (igraph::is_connected(G)) break
  }
  G
}

sbm_pm <- function(ncomm = 2) {
  p_within  <- runif(ncomm, 0.5, 0.9)
  p_between <- runif(ncomm, 0.1, 0.4)
  mat <- matrix(0, ncomm, ncomm); diag(mat) <- p_within
  mat[lower.tri(mat)] <- runif(1, min(p_between), max(p_between))
  mat + t(mat) - diag(diag(mat))
}

gen <- list(
  ER = function(n) { p <- sample(ER_P, 1)
  list(g = igraph::sample_gnp(n, p), param = paste0("p=", p)) },
  Spatial = function(n) { r <- sample(SP_R, 1)
  list(g = fastSpatialNetwork(n, r, makeConnected = TRUE), param = paste0("r=", r)) },
  SF = function(n) { pw <- sample(SF_POWER, 1); m <- sample(SF_M[SF_M <= n - 1], 1)
  list(g = igraph::sample_pa(n, power = pw, m = m, directed = FALSE, algorithm = "psumtree"),
       param = paste0("power=", pw, ",m=", m)) },
  SW = function(n) { p <- sample(SW_P, 1); nei <- sample(SW_NEI, 1)
  list(g = igraph::sample_smallworld(dim = 1, size = n, nei = nei, p = p),
       param = paste0("p=", p, ",nei=", nei)) },
  sbm = function(n) { b1 <- round(0.4 * n)
  list(g = igraph::sample_sbm(n, pref.matrix = sbm_pm(2), block.sizes = c(b1, n - b1)),
       param = "2 blocks") }
)

draw_connected <- function(f, n) {
  for (t in seq_len(MAX_TRIES)) {
    d <- f(n); g <- igraph::simplify(d$g)
    if (igraph::ecount(g) > 0 && igraph::is_connected(g)) return(list(g = g, param = d$param, tries = t))
  }
  list(g = NULL, param = NA_character_, tries = MAX_TRIES)
}

set.seed(2026)   # same seed as the 2026-10-10 runs: results reproduce exactly
rows <- list(); draws_log <- list()
for (n in SIZES) {
  for (cl in names(gen)) {
    draws <- lapply(seq_len(N_PER), function(i) draw_connected(gen[[cl]], n))
    ok    <- !vapply(draws, function(d) is.null(d$g), logical(1))
    draws_log[[length(draws_log) + 1]] <- tibble::tibble(
      true = cl, n = n, n_kept = sum(ok), n_gave_up = sum(!ok),
      mean_tries = round(mean(vapply(draws, `[[`, numeric(1), "tries")), 2))
    if (!any(ok)) next
    graphs <- lapply(draws[ok], `[[`, "g")
    
    gf <- calc_graph_features(graphs)
    gf$graph_name <- paste(cl, n, seq_along(graphs), sep = "_")
    fc <- gf %>% dplyr::select(dplyr::any_of(FEATURE_COLS)) %>%
      dplyr::mutate_if(is.character, factor) %>% janitor::clean_names() %>% tibble::as_tibble()
    
    rows[[length(rows) + 1]] <- dplyr::bind_cols(
      tibble::tibble(true = cl, n = n, param = vapply(draws[ok], `[[`, "", "param"),
                     density = 2 * gf$edges / (gf$order * (gf$order - 1))),
      summarise_probs(prob_main(model, fc)))
    message(sprintf("n = %3d  %-7s  kept %d", n, cl, sum(ok)))
  }
}
val <- dplyr::bind_rows(rows) %>%
  dplyr::mutate(true = factor(true, CLASS_ORDER), predicted = factor(predicted, CLASS_ORDER),
                correct = true == predicted)
draws_log <- dplyr::bind_rows(draws_log)
cm <- with(val, table(true, predicted))

cat("\nOverall accuracy:", round(100 * mean(val$correct), 1), "% of", nrow(val), "networks\n")
cat("\nConfusion matrix, rows = true:\n"); print(cm)
cat("\nWhat each true class is called (row %):\n"); print(round(100 * prop.table(cm, 1), 1))

rp <- tibble::tibble(class     = CLASS_ORDER,
                     recall    = round(100 * diag(cm) / rowSums(cm), 1),
                     precision = round(100 * diag(cm) / colSums(cm), 1))
cat("\nRecall and precision (precision assumes equal class numbers):\n"); print(rp)

acc_size <- val %>% dplyr::group_by(true, n) %>%
  dplyr::summarise(acc = round(mean(correct), 3), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = n, values_from = acc)
cat("\nAccuracy by true class x size:\n"); print(acc_size, width = Inf)

acc_density <- val %>%
  dplyr::mutate(density_band = cut(density, c(0, .1, .25, .5, .75, 1), include.lowest = TRUE)) %>%
  dplyr::group_by(density_band) %>%
  dplyr::summarise(accuracy = round(mean(correct), 3), n_graphs = dplyr::n())
cat("\nAccuracy by density band:\n"); print(acc_density)

conf_val <- val %>% dplyr::group_by(correct) %>%
  dplyr::summarise(median_p_max = round(stats::median(p_max), 3),
                   median_margin = round(stats::median(margin), 3), n_graphs = dplyr::n())
cat("\nConfidence of correct vs incorrect calls:\n"); print(conf_val)

cat("\nRedraws needed to get a connected graph:\n"); print(draws_log, n = Inf)

utils::write.csv(cbind(true = rownames(cm), as.data.frame.matrix(cm), rp[, c("recall", "precision")]),
                 file.path(OUT_DIR, "tableS_point2b_confusion.csv"), row.names = FALSE)
utils::write.csv(acc_size,    file.path(OUT_DIR, "tableS_point2b_accuracy_by_size.csv"), row.names = FALSE)
utils::write.csv(acc_density, file.path(OUT_DIR, "tableS_point2b_accuracy_by_density.csv"), row.names = FALSE)
utils::write.csv(val,         file.path(OUT_DIR, "point2_validation_small_networks.csv"), row.names = FALSE)
utils::write.csv(draws_log,   file.path(OUT_DIR, "point2_validation_draws.csv"), row.names = FALSE)


# =============================================================================
# (e) EMPIRICAL CLASSIFICATION UNCERTAINTY (same inputs as script 02)
# =============================================================================
hr("(e) Classification uncertainty for the empirical networks")

emp <- utils::read.csv(here::here("data", "GraphFeatOnAllAnimNets.csv")) %>%
  dplyr::rename(graph_name = GraphNames) %>%
  dplyr::select(dplyr::all_of(FEATURE_COLS)) %>%
  janitor::clean_names() %>%
  dplyr::filter(order > 10) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.integer), as.double)) %>%
  tibble::as_tibble()

emp_prepped <- recipes::recipe(graph_name ~ ., data = emp) %>%
  recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
  recipes::prep() %>% recipes::juice()

P_emp   <- prob_main(model, emp_prepped)
emp_unc <- dplyr::bind_cols(
  tibble::tibble(graph_name = as.character(emp_prepped$graph_name),
                 order = emp_prepped$order,
                 size_band = cut(emp_prepped$order, c(10, 49, 1000, Inf),
                                 labels = c("<50", "50-1000", ">1000"))),
  summarise_probs(P_emp), tibble::as_tibble(round(P_emp, 4)))

cat("Predicted classes (expect Spatial 392, SF 167, ER 30, SW 24, sbm 23):\n")
print(table(emp_unc$predicted))

unc_band <- emp_unc %>% dplyr::group_by(size_band) %>%
  dplyr::summarise(n = dplyr::n(), median_p_max = round(stats::median(p_max), 3),
                   pct_p_max_lt_0.5 = round(100 * mean(p_max < 0.5), 1),
                   median_margin = round(stats::median(margin), 3))
unc_class <- emp_unc %>% dplyr::group_by(predicted) %>%
  dplyr::summarise(n = dplyr::n(), n_below_50 = sum(order < 50),
                   median_p_max = round(stats::median(p_max), 3),
                   pct_p_max_lt_0.5 = round(100 * mean(p_max < 0.5), 1))
cat("\nBy size band:\n"); print(unc_band)
cat("\nBy predicted class:\n"); print(unc_class)

utils::write.csv(emp_unc, file.path(OUT_DIR, "dataS_point2e_empirical_class_probabilities.csv"), row.names = FALSE)
utils::write.csv(dplyr::bind_rows(dplyr::mutate(unc_band,  group = as.character(size_band)),
                                  dplyr::mutate(unc_class, group = predicted)) %>%
                   dplyr::select(group, dplyr::everything(), -dplyr::any_of(c("size_band", "predicted"))),
                 file.path(OUT_DIR, "tableS_point2e_uncertainty_summary.csv"), row.names = FALSE)


# =============================================================================
# SUMMARY: the answer to each part of point 2, from the numbers above
# =============================================================================
hr("SUMMARY")

tr_all <- dplyr::filter(size_tab, class == "All")
cat(sprintf(
  "(a) Training networks: %s, order %g-%g; networks below 50 nodes: %d.
    Empirical networks below 50 nodes: %d of %d -> outside the training size range.\n",
  format(tr_all$n_networks, big.mark = ","), tr_all$min_order, tr_all$max_order,
  tr_all$n_order_lt_50, n_emp_lt50, nrow(emp_fig)))

cat(sprintf(
  "\n(b) Existing classifier on %d simulated networks, n = 11-45: overall accuracy %.1f%%.
    Recall: %s.
    Accuracy by size and density: tableS_point2b_accuracy_by_size.csv, _by_density.csv.\n",
  nrow(val), 100 * mean(val$correct),
  paste(sprintf("%s %.1f%%", rp$class, rp$recall), collapse = ", ")))

top_out <- utils::head(dplyr::filter(range_tab, n_empirical_outside > 0), 3)
cat(sprintf(
  "\n(c) Inside the training range on all %d features: %d of %d networks
    (<50 nodes: %d of %d; 50-1000 nodes: %d of %d).
    Most often outside: %s.\n",
  length(feat_names), sum(all_in), nrow(emp_fig),
  all_tab$n_inside_all_features[all_tab$group == "<50"],     all_tab$n[all_tab$group == "<50"],
  all_tab$n_inside_all_features[all_tab$group == "50-1000"], all_tab$n[all_tab$group == "50-1000"],
  paste(sprintf("%s (%d)", top_out$feature, top_out$n_empirical_outside), collapse = ", ")))

cat(sprintf(
  "\n(d) %.1f%% = %d of %d empirical networks inside the convex hull of the %s
    training networks in the PC1-PC2 plane (PC1 %.1f%% + PC2 %.1f%% = %.1f%% of variance).\n",
  cov_tab$pct_inside_convex_hull, cov_tab$n_inside_convex_hull, cov_tab$n_empirical,
  format(cov_tab$n_simulated, big.mark = ","), pc1_var, pc2_var, pc1_var + pc2_var))

cat(sprintf(
  "\n(e) Median top-class probability: %s.
    In the simulation (b), wrong calls had median top-class probability %.3f
    (correct calls %.3f), so below 50 nodes the validation error rates are the
    uncertainty measure. Per-network probabilities: dataS_point2e_empirical_class_probabilities.csv.\n",
  paste(sprintf("%s nodes %.3f (n = %d)", unc_band$size_band, unc_band$median_p_max, unc_band$n), collapse = ", "),
  conf_val$median_p_max[!conf_val$correct], conf_val$median_p_max[conf_val$correct]))

message("PCA_Applicability.R complete.")