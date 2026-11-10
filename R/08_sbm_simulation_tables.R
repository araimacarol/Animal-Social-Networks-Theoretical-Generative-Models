# =============================================================================
# 08_sbm_simulation_tables.R
# -----------------------------------------------------------------------------
#Supplementary Tables S16 and S17 (SBM simulations with
# lower / higher between-community connection probability).
#
# The simulations were run with generate_sbm_networks() in
# sbm-sim-different-communities.R, which wrote one row per simulated network
# (features and predicted class from the same network) to:
#   outputs/all_sbm_preds_lower_between_connectn_probs.csv   -> Table S16
#   outputs/all_sbm_preds_higher_between_connectn_probs.csv  -> Table S17
# This script writes the LaTeX
# tables directly from the CSV files, after checking each row
# (mean degree == 2E/N), and summarises the predicted classes.
#
# Design (from generate_sbm_networks()): k = 2-6 communities, n = 100-1000 in
# steps of 50, one network per k and n; block sizes floor(n/k) with the
# remainder added to block 1; within-community p ~ U(0.5, 0.8) per community;
# between-community q drawn per pair of communities from U(0.01, 0.04) (lower)
# or U(0.1, 0.4) (higher). Networks classified with the first-stage classifier.
#
# Outputs: outputs/tableS16_sbm_lower.tex, outputs/tableS17_sbm_higher.tex,
#          outputs/tableS_point3_sbm_summary.csv
#
#   source(here::here("R", "08_sbm_simulation_tables.R"))
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

FILES <- c(
  lower  = here::here("outputs", "all_sbm_preds_lower_between_connectn_probs.csv"),
  higher = here::here("outputs", "all_sbm_preds_higher_between_connectn_probs.csv"))

sims <- dplyr::bind_rows(lapply(names(FILES), function(sc) {
  d <- utils::read.csv(FILES[[sc]], check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[names(d) == ".pred_class"] <- "predicted_class"
  names(d)[names(d) == "k_communities"] <- "k"
  d$scenario <- sc
  d
}))

## ---- Row checks -------------------------------------------------------------
stopifnot(all(abs(sims$mean_degree - 2 * sims$edges / sims$order) < 1e-6))
stopifnot(all(sims$order == sims$n_nodes))
cat("Networks per scenario:\n"); print(table(sims$scenario))

## ---- The stored predictions come from the recovered first-stage model ------
## Re-classify every stored feature row with the recovered booster and check
## that it gives the class written in the CSV.
feat_cols <- janitor::make_clean_names(FEATURE_COLS)
re_pred <- predict_main_model(load_main_model(),
                              tibble::as_tibble(sims[, feat_cols]) %>%
                                dplyr::mutate(graph_name = factor(graph_name)))$predicted_class
cat("Stored predictions reproduced by the recovered model:",
    sum(re_pred == sims$predicted_class), "of", nrow(sims), "\n")
stopifnot(all(re_pred == sims$predicted_class))

## ---- Summary: predicted classes and modularity -----------------------------
summ <- sims %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(
    networks        = dplyr::n(),
    mean_modularity = round(mean(modularity), 3),
    n_sbm           = sum(predicted_class == "sbm"),
    n_geometric     = sum(predicted_class == "Spatial"),
    n_small_world   = sum(predicted_class == "SW"),
    n_other         = sum(!predicted_class %in% c("sbm", "Spatial", "SW")),
    .groups = "drop")
print(as.data.frame(summ), row.names = FALSE)
utils::write.csv(summ, here::here("outputs", "tableS_point3_sbm_summary.csv"), row.names = FALSE)

## ---- LaTeX tables -------------------------------------------------------------
fmt <- function(x, d) formatC(x, format = "f", digits = d)
write_sbm_table <- function(d, caption, label, file) {
  d <- d %>% dplyr::arrange(order, k)
  hdr <- paste("k & Pred & Order & Edges & M.Ecc & M.Path & Energy & Mod & Diam & Betw & Trans &",
               "Spec.R & Eig.C & Deg.C & M.Deg & MinCut & Fied & N.Fied & Close & Assort \\\\")
  body <- sprintf(
    "%d & %s & %d & %d & %s & %s & %s & %s & %d & %s & %s & %s & %s & %s & %s & %d & %s & %s & %s & %s \\\\",
    as.integer(d$k), d$predicted_class, as.integer(d$order), as.integer(d$edges),
    fmt(d$mean_eccentr, 3), fmt(d$mean_path_length, 3), fmt(d$graph_energy, 1),
    fmt(d$modularity, 3), as.integer(d$diameter), fmt(d$betw_centr, 4), fmt(d$transitivity, 3),
    fmt(d$spectral_radius, 2), fmt(d$eigen_centr, 3), fmt(d$deg_centr, 3), fmt(d$mean_degree, 2),
    as.integer(d$min_cut), fmt(d$fiedler_value, 2), fmt(d$normalized_fiedler_value, 3),
    fmt(d$closeness_centr, 5), fmt(d$deg_assort_coef, 3))
  head <- c("\\toprule", hdr, "\\midrule")
  lines <- c(
    "{\\tiny", "\\setlength{\\tabcolsep}{1.7pt}",
    "\\begin{longtable}{@{}rlrrrrrrrrrrrrrrrrrr@{}}",
    paste0("\\caption{\\changeto{}{", caption, "}}"), paste0("\\label{", label, "}\\\\"),
    head, "\\endfirsthead", head, "\\endhead", "\\bottomrule", "\\endfoot",
    body, "\\end{longtable}", "}")
  writeLines(lines, file); message("Written: ", file)
}

caption_text <- function(range_txt, s) {
  paste0("Simulated stochastic block models with ", range_txt,
         " between-community connection probability, one network per number of communities ($k$) and size. ",
         "Features and predicted class are computed from the same network. ",
         "Predicted classes: SBM ", s$n_sbm, ", geometric graph ", s$n_geometric,
         ", small-world ", s$n_small_world, " of ", s$networks,
         " networks; mean modularity ", fmt(s$mean_modularity, 3), ". ",
         "Pred = predicted class (sbm = stochastic block model, Spatial = geometric graph, SW = small-world); ",
         "M.Ecc = mean eccentricity; M.Path = mean path length; Energy = graph energy; Mod = modularity; ",
         "Diam = diameter; Betw = betweenness centralisation; Trans = transitivity; Spec.R = spectral radius; ",
         "Eig.C = eigenvector centralisation; Deg.C = degree centralisation; M.Deg = mean degree ($2E/N$); ",
         "MinCut = minimum cut; Fied = Fiedler value; N.Fied = normalised Fiedler value; ",
         "Close = mean closeness centrality; Assort = degree assortativity.")
}
s_low  <- dplyr::filter(summ, scenario == "lower")
s_high <- dplyr::filter(summ, scenario == "higher")
write_sbm_table(dplyr::filter(sims, scenario == "lower"),
  caption_text("lower ($0.01 \\leq q \\leq 0.04$)", s_low),
  "tab:sbm_summary", here::here("outputs", "tableS16_sbm_lower.tex"))
write_sbm_table(dplyr::filter(sims, scenario == "higher"),
  caption_text("higher ($0.1 \\leq q \\leq 0.4$)", s_high),
  "tab:sbm_summary2", here::here("outputs", "tableS17_sbm_higher.tex"))

message("08_sbm_simulation_tables.R complete.")
