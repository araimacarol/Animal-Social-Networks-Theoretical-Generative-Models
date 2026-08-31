# =============================================================================
# 01_functions.R  --  all reusable functions for the pipeline
# -----------------------------------------------------------------------------
# Pure functions only: nothing here reads/writes files at source() time.
# The numbered pipeline scripts call these. Grouped into:
#   A. Graph feature extraction
#   B. Recovered main (first-stage) model: load + predict + classify a graph
#   C. Group-level module classification
#   D. Community-detection comparison across algorithms
#   E. Grouped (leave-one-group-out) cross-validation for the metadata models
#   F. Metadata-model recovery + evaluation + LaTeX tables
#   G. Threshold-sensitivity test for unweighted .edges networks
#   H. Plot helpers (metadata bars, SHAP)
# =============================================================================


# =============================================================================
# 0. EDGE-FILE RESOLUTION
# -----------------------------------------------------------------------------
# Edge lists keep their original download names; some end ".csv .edges" (with a
# space). Callers pass a bare network name and these helpers find the real file.
# =============================================================================

resolve_edge_file <- function(name, dir = dir_edges) {
  candidates <- file.path(dir, paste0(name, c(".edges", ".csv .edges", ".csv.edges")))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

resolve_edge_files <- function(names, dir = dir_edges) {
  paths   <- vapply(names, resolve_edge_file, character(1), dir = dir)
  missing <- names[is.na(paths)]
  if (length(missing) > 0) {
    stop("No edge file found in ", dir, " for: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  unname(paths)
}


# =============================================================================
# A. GRAPH FEATURE EXTRACTION
# =============================================================================

#' Normalised Laplacian of a graph
normalized_laplacian <- function(g) igraph::laplacian_matrix(g,
                            normalization = c("symmetric"))

#' Compute the topological feature vector used by the classifier.
#'
#' @param graphs A list of igraph objects.
#' @return A data frame, one row per graph, with the feature columns the
#'   first-stage classifier was trained on.
calc_graph_features <- function(graphs = NULL) {
  base_features <- c(
    "order", "edges", "connected", "max_component", "minDegree", "maxDegree",
    "mean_degree", "minCut", "FiedlerValue", "Normalized_FiedlerValue",
    "closeness_centr", "modularity", "diameter", "betw_centr", "transitivity",
    "threshold", "spectral_radius"
  )
  df <- as.data.frame(matrix(ncol = length(base_features), nrow = length(graphs)))
  colnames(df) <- base_features
  
  ## Vectorizable, no objects needed
  df$order        <- as.numeric(lapply(graphs, igraph::gorder))
  df$edges        <- as.numeric(lapply(graphs, igraph::gsize))
  df$connected    <- as.numeric(lapply(graphs, igraph::is_connected))
  df$minCut       <- as.numeric(lapply(graphs, igraph::min_cut))
  df$diameter     <- as.numeric(lapply(graphs, igraph::diameter))
  df$transitivity <- as.numeric(lapply(graphs, igraph::transitivity))
  
  deg            <- lapply(graphs, igraph::degree)
  df$minDegree   <- as.numeric(lapply(deg, min))
  df$maxDegree   <- as.numeric(lapply(deg, max))
  df$mean_degree <- as.numeric(lapply(deg, mean))
  
  ## Objects, filled per-graph in a loop
  communities <- lapply(graphs, igraph::cluster_walktrap)
  Adj         <- lapply(graphs, igraph::as_adjacency_matrix, sparse = TRUE)
  L           <- lapply(graphs, igraph::laplacian_matrix)
  Norm_Lap    <- lapply(graphs, normalized_laplacian)
  
  for (i in seq_along(graphs)) {
    if (is.null(graphs[[i]]$type)) graphs[[i]]$type <- "untyped"
    
    df$modularity[i]      <- igraph::modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric = TRUE, only.values = TRUE)$values[1]
    
    ev_L  <- eigen(L[[i]],        symmetric = TRUE, only.values = TRUE)$values
    ev_NL <- eigen(Norm_Lap[[i]], symmetric = TRUE, only.values = TRUE)$values
    df$FiedlerValue[i]            <- ev_L[length(ev_L) - 1]
    df$Normalized_FiedlerValue[i] <- ev_NL[length(ev_NL) - 1]
    
    df$eigen_centr[i]        <- igraph::centr_eigen(graphs[[i]])$centralization
    df$deg_centr[i]          <- igraph::centr_degree(graphs[[i]])$centralization
    df$betw_centr[i]         <- igraph::centr_betw(graphs[[i]])$centralization
    df$max_component[i]      <- max(igraph::components(graphs[[i]])$csize)
    df$mean_eccentr[i]       <- mean(igraph::eccentricity(graphs[[i]]))
    df$radius[i]             <- igraph::radius(graphs[[i]])
    df$mean_path_length[i]   <- igraph::mean_distance(graphs[[i]])
    df$graph_energy[i]       <- sum(abs(eigen(Adj[[i]], symmetric = TRUE, only.values = TRUE)$values))
    df$deg_assort_coef[i]    <- igraph::assortativity_degree(graphs[[i]])
    df$threshold[i]          <- 1 / df$spectral_radius[i]
    
    df$closeness_centr[i] <- if (isTRUE(df$connected[i] == 1)) {
      mean(igraph::closeness(graphs[[i]]))
    } else {
      -1  # sentinel: not defined for disconnected graphs
    }
  }
  df
}

## The 19 predictor columns (post janitor::clean_names) the model expects.
FEATURE_COLS <- c(
  "graph_name", "order", "edges", "mean_eccentr", "mean_path_length",
  "graph_energy", "modularity", "diameter", "betw_centr", "transitivity",
  "spectral_radius", "eigen_centr", "deg_centr", "mean_degree", "minCut",
  "FiedlerValue", "Normalized_FiedlerValue", "closeness_centr", "deg_assort_coef"
)


# =============================================================================
# B. RECOVERED FIRST-STAGE MODEL: load, predict, classify one graph
# -----------------------------------------------------------------------------
# The original workflow object saved from tidymodels became unloadable after an
# xgboost version bump ("object is corrupted or is from an incompatible XGBoost
# version"). We therefore ship a version-portable JSON booster + the trained
# recipe, and predict with the native booster. See README > "Model recovery".
# =============================================================================

#' Load the recovered first-stage model (portable booster + trained recipe).
load_main_model <- function(booster_json = here::here("models", "booster_portable_mainmodel.json"),
                            recipe_rds   = here::here("models", "preproc_trained_mainmodel.rds"),
                            class_order  = CLASS_ORDER) {
  if (!file.exists(booster_json)) stop("Recovered booster not found: ", booster_json)
  if (!file.exists(recipe_rds))   stop("Trained recipe not found: ",   recipe_rds)
  list(
    booster     = xgboost::xgb.load(booster_json),
    preproc     = readRDS(recipe_rds),
    class_order = class_order
  )
}

#' Predict class labels from a feature data frame with the recovered booster.
predict_main_model <- function(model, feature_df) {
  feature_df <- feature_df %>% dplyr::mutate(dplyr::across(where(is.integer), as.double))
  baked <- recipes::bake(model$preproc, new_data = feature_df)
  
  drop_cols     <- intersect(c("graph_name", "target", "class", ".pred_class"), names(baked))
  correct_order <- setdiff(names(baked), drop_cols)
  X <- as.matrix(baked[, correct_order, drop = FALSE]); storage.mode(X) <- "double"
  
  # xgboost >= 2.0 dropped `reshape`; multi:softprob already returns an n x K matrix.
  raw <- stats::predict(model$booster, X)
  if (!is.matrix(raw)) raw <- matrix(raw, nrow = nrow(X), byrow = TRUE)  # older xgboost
  idx <- if (is.matrix(raw) && ncol(raw) > 1) {
    max.col(raw, ties.method = "first")            # multi:softprob
  } else {
    as.integer(round(as.numeric(raw))) + 1L        # multi:softmax (0-based)
  }
  tibble::tibble(predicted_class = model$class_order[idx])
}

#' Compute features for one graph and classify it (returns a class string or NULL).
classify_graph <- function(g, model, tag = "graph") {
  if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) return(NULL)
  feats <- calc_graph_features(list(g)); feats$graph_name <- tag
  feats_clean <- feats %>%
    dplyr::filter(connected == 1) %>%
    dplyr::select(dplyr::any_of(FEATURE_COLS)) %>%
    dplyr::mutate_if(is.character, factor) %>%
    janitor::clean_names() %>%
    tibble::as_tibble()
  if (nrow(feats_clean) == 0) return(NULL)
  pred <- tryCatch(predict_main_model(model, feats_clean)$predicted_class,
                   error = function(e) { message("  predict failed: ", conditionMessage(e)); NA_character_ })
  if (length(pred) == 0 || is.na(pred[1])) return(NULL)
  as.character(pred[1])
}


# =============================================================================
# C. GROUP-LEVEL MODULE CLASSIFICATION
# =============================================================================

#' Split one network into its communities and classify each community.
#'
#' @param file Path to an .edges file.
#' @param badger Logical; if TRUE use the prescribed badger group memberships
#'   from Badger_Groups.csv, otherwise detect communities with walktrap.
#' @param model A model object from load_main_model().
#' @param badger_groups_file Path to the prescribed-groups CSV.
group_modules <- function(file, badger = FALSE, model = load_main_model(),
                          badger_groups_file = here::here("data", "Badger_Groups.csv")) {
  df <- utils::read.table(file, sep = " ", fill = TRUE, header = FALSE)
  g  <- igraph::graph_from_data_frame(as.matrix(df), directed = FALSE)
  g  <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g  <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
  
  comps <- igraph::components(g, mode = "weak")
  graph <- igraph::induced_subgraph(g, igraph::V(g)[comps$membership == which.max(comps$csize)])
  
  if (badger) {
    bg <- utils::read.csv(badger_groups_file)
    igraph::V(graph)$membership <- bg$Group[match(igraph::V(graph)$name, bg$NodeID)]
    membership <- igraph::V(graph)$membership
    mod_val    <- igraph::modularity(graph, igraph::membership(igraph::cluster_walktrap(graph)))
  } else {
    cw         <- igraph::cluster_walktrap(graph)
    membership <- cw$membership
    mod_val    <- igraph::modularity(graph, membership)
  }
  
  ## connected subgraphs, keeping original group IDs as names
  subgraphs <- list()
  for (grp in unique(membership)) {
    nodes <- which(membership == grp)
    sub   <- igraph::induced_subgraph(graph, nodes)
    if (igraph::is_connected(sub) && igraph::ecount(sub) > 0)
      subgraphs[[as.character(grp)]] <- sub
  }
  
  module_summary <- vapply(seq_along(subgraphs), function(i) {
    sprintf("Group %s has %d nodes and %d edges",
            names(subgraphs)[i],
            as.integer(igraph::vcount(subgraphs[[i]])),
            as.integer(igraph::ecount(subgraphs[[i]])))
  }, character(1))
  
  features_df <- do.call(rbind, lapply(subgraphs, function(s) calc_graph_features(list(s))))
  features_df <- features_df %>% dplyr::filter(connected == 1)
  features_df$graph_name <- names(subgraphs)[seq_len(nrow(features_df))]
  
  features_clean <- features_df %>%
    dplyr::select(dplyr::any_of(FEATURE_COLS)) %>%
    dplyr::mutate_if(is.character, factor) %>%
    janitor::clean_names() %>%
    tibble::as_tibble()
  
  predictions <- predict_main_model(model, features_clean) %>%
    dplyr::bind_cols(features_clean)
  
  list(
    Modules                 = module_summary,
    subgraphs_igraph        = subgraphs,
    Graph_features          = features_df,
    Class_Predictions       = predictions,
    Complete_connected_graph = graph,
    Modularity_Val          = mod_val
  )
}


# =============================================================================
# D. COMMUNITY-DETECTION COMPARISON ACROSS ALGORITHMS
# -----------------------------------------------------------------------------
# Are group-level results stable across community-detection methods?
# We run several algorithms + the prescribed grouping, classifies every
# community >= MIN_GROUP_SIZE, and records per-network group sizes and the
# predicted generative-model class. Produces merged_community_predictions.csv.
# =============================================================================

#' For each network and method, report group sizes and the predicted class of
#' every community that clears `min_group_size`.
summarize_group_sizes_and_predictions <- function(
    edge_files,
    model              = load_main_model(),
    methods            = c("walktrap", "louvain", "leading_eigen", "infomap", "prescribed"),
    badger_groups_file = here::here("data", "Badger_Groups.csv"),
    min_group_size     = MIN_GROUP_SIZE) {
  
  badger_groups <- if (file.exists(badger_groups_file)) utils::read.csv(badger_groups_file) else NULL
  edge_files    <- edge_files[file.exists(edge_files)]
  
  build_base_graph <- function(file) {
    df <- utils::read.table(file, sep = " ", fill = TRUE, header = FALSE)
    g  <- igraph::graph_from_data_frame(as.matrix(df), directed = FALSE)
    g  <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    g  <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    comps <- igraph::components(g, mode = "weak")
    igraph::induced_subgraph(g, igraph::V(g)[comps$membership == which.max(comps$csize)])
  }
  
  ## Records (graph, method) pairs where the algorithm failed to produce a valid
  ## partition -- currently only leading-eigenvector / ARPACK non-convergence.
  failures <<- tibble::tibble()
  
  get_membership <- function(graph, method, is_badger, graph_name = NA_character_) {
    if (method == "prescribed") {
      if (is_badger && !is.null(badger_groups))
        return(badger_groups$Group[match(igraph::V(graph)$name, badger_groups$NodeID)])
      return(NULL)
    }
    
    converged <- TRUE
    comm <- withCallingHandlers(
      tryCatch(switch(method,
                      walktrap      = igraph::cluster_walktrap(graph),
                      louvain       = igraph::cluster_louvain(graph),
                      leading_eigen = igraph::cluster_leading_eigen(graph),
                      infomap       = igraph::cluster_infomap(graph),
                      stop("Unknown method: ", method)), error = function(e) NULL),
      warning = function(w) {
        ## ARPACK bails out but igraph still RETURNS an object. That partition is
        ## not a valid result and must not be classified.
        if (grepl("ARPACK|converge", conditionMessage(w), ignore.case = TRUE)) {
          converged <<- FALSE
          invokeRestart("muffleWarning")
        }
      })
    
    if (is.null(comm) || !converged) {
      message("  ", method, " failed on ", graph_name,
              if (!converged) " (ARPACK did not converge)" else " (error)", " -- excluded")
      failures <<- dplyr::bind_rows(failures, tibble::tibble(
        graph_name = graph_name, method = method,
        reason = if (!converged) "arpack_no_convergence" else "algorithm_error"))
      return(NULL)
    }
    igraph::membership(comm)
  }
  
  extract_subgraphs <- function(graph, membership_vec, min_size) {
    out <- list()
    for (grp in unique(membership_vec)) {
      nodes <- which(membership_vec == grp)
      if (length(nodes) < min_size) next
      sub <- igraph::induced_subgraph(graph, nodes)
      if (igraph::is_connected(sub) && igraph::ecount(sub) > 0) out[[as.character(grp)]] <- sub
    }
    out
  }
  
  classify_subgraphs <- function(subgraph_list, graph_name, method) {
    if (length(subgraph_list) == 0) return(tibble::tibble())
    feats <- calc_graph_features(subgraph_list)
    feats$group_id <- names(subgraph_list)
    # igraph >= 2.0: vcount() returns a double, so the vapply template must be
    # numeric(1); cast back to integer for the group-size threshold + output.
    feats$n_nodes  <- as.integer(vapply(subgraph_list, igraph::vcount, numeric(1)))
    feats_clean <- feats %>%
      dplyr::filter(connected == 1) %>%
      dplyr::rename(graph_name_col = group_id) %>%
      dplyr::select(dplyr::any_of(c(FEATURE_COLS[-1], "graph_name_col", "n_nodes"))) %>%
      dplyr::mutate_if(is.character, factor) %>%
      janitor::clean_names() %>%
      dplyr::mutate(dplyr::across(where(is.integer), as.double)) %>%
      tibble::as_tibble()
    if (nrow(feats_clean) == 0) return(tibble::tibble())
    pred <- predict_main_model(model, feats_clean)$predicted_class
    tibble::tibble(graph_name = graph_name, method = method,
                   group_id = feats_clean$graph_name_col,
                   n_nodes  = feats_clean$n_nodes, predicted_class = pred)
  }
  
  size_summary <- tibble::tibble(); prediction_details <- tibble::tibble()
  for (file in edge_files) {
    graph_name <- tools::file_path_sans_ext(basename(file))
    is_badger  <- grepl("badger", graph_name, ignore.case = TRUE)
    g <- tryCatch(build_base_graph(file), error = function(e) NULL)
    if (is.null(g)) next
    total_nodes <- as.integer(igraph::vcount(g))
    
    for (method in methods) {
      mem <- get_membership(g, method, is_badger, graph_name)
      if (is.null(mem)) next
      sizes <- as.integer(table(mem))
      size_summary <- dplyr::bind_rows(size_summary, tibble::tibble(
        graph_name = graph_name, method = method,
        total_nodes_in_network = total_nodes,
        n_groups_detected = length(sizes),
        group_sizes = paste(sort(sizes, decreasing = TRUE), collapse = ", "),
        n_groups_passing_threshold = sum(sizes >= min_group_size),
        pct_groups_classifiable = round(100 * sum(sizes >= min_group_size) / length(sizes), 1)))
      subgraphs <- extract_subgraphs(g, mem, min_group_size)
      prediction_details <- dplyr::bind_rows(
        prediction_details, classify_subgraphs(subgraphs, graph_name, method))
    }
  }
  list(size_summary = size_summary, predictions = prediction_details,
       failures = failures)
}

#' Join size summary + predictions into the merged table
merge_community_predictions <- function(cd) {
  cd$size_summary %>%
    dplyr::left_join(
      cd$predictions %>% dplyr::rename(n_nodes_from_passed_group = n_nodes),
      by = c("graph_name", "method")) %>%
    dplyr::select(graph_name, method, total_nodes_in_network, n_groups_detected,
                  group_sizes, n_groups_passing_threshold, group_id,
                  n_nodes_from_passed_group, predicted_class) %>%
    dplyr::arrange(graph_name, method, group_id) %>%
    dplyr::mutate(predicted_class_label = dplyr::recode(predicted_class,
                                                        "ER" = "Erd\\H{o}s-R\\'enyi", "sbm" = "SBM", "SF" = "Scale Free",
                                                        "Spatial" = "Geometric Graph", "SW" = "Small World"))
}


# =============================================================================
# E. GROUPED (LEAVE-ONE-GROUP-OUT) CROSS-VALIDATION
# -----------------------------------------------------------------------------
# Does a random 70:30 split overestimates performance because many
# networks share a species / study design / protocol. This section holds out an entire group
# at a time and re-evaluates. Produces pipeline_overall_performance2.csv and
# pipeline_classwise_performance2.csv.
# =============================================================================

run_ecological_network_pipeline <- function(
    data_path = here::here("models", "df_final_meta_data.rds"),
    seed = 123, top_n_species = 5, tuning_grid_size = 15) {
  
  df_raw <- readRDS(data_path) %>%
    dplyr::mutate(dplyr::across(c(target, species, class, interaction_type), as.factor))
  set.seed(seed)
  df_shuffled <- df_raw[sample(seq_len(nrow(df_raw))), ]
  
  xgb_spec <- parsnip::boost_tree(
    mtry = tune::tune(), trees = tune::tune(),
    learn_rate = tune::tune(), tree_depth = tune::tune()) %>%
    parsnip::set_engine("xgboost") %>%
    parsnip::set_mode("classification")
  
  overall <- tibble::tibble(); classwise <- tibble::tibble()
  
  ## overall + class-wise + majority-baseline metrics from one confusion matrix
  extract_matrix_metrics <- function(cm, model_name, validation_type, held_out) {
    class_names <- rownames(cm)
    cw <- purrr::map_df(class_names, function(cl) {
      tp <- cm[cl, cl]; tot_pred <- sum(cm[, cl]); tot_act <- sum(cm[cl, ])
      precision <- if (tot_pred == 0) 0 else tp / tot_pred
      recall    <- if (tot_act  == 0) 0 else tp / tot_act
      f1        <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
      tibble::tibble(Model = model_name, Validation = validation_type, Held_Out_Group = held_out,
                     Target_Class = cl, Class_Precision = round(precision, 4),
                     Class_Recall = round(recall, 4), Class_F1 = round(f1, 4),
                     Class_Support = tot_act)
    })
    total <- sum(cm); w <- cw$Class_Support / sum(cw$Class_Support)
    majority <- class_names[which.max(rowSums(cm))]
    ov <- tibble::tibble(
      Model = model_name, Validation = validation_type, Held_Out_Group = held_out,
      Overall_Accuracy   = round(sum(diag(cm)) / total, 4),
      Balanced_Accuracy  = round(mean(cw$Class_Recall, na.rm = TRUE), 4),
      Macro_Precision    = round(mean(cw$Class_Precision, na.rm = TRUE), 4),
      Macro_F1           = round(mean(cw$Class_F1, na.rm = TRUE), 4),
      Weighted_Precision = round(sum(cw$Class_Precision * w, na.rm = TRUE), 4),
      Weighted_Recall    = round(sum(cw$Class_Recall    * w, na.rm = TRUE), 4),
      Weighted_F1        = round(sum(cw$Class_F1         * w, na.rm = TRUE), 4),
      Baseline_Accuracy  = round(sum(cm[majority, ]) / total, 4),
      Total_Test_Samples = total)
    list(overall = ov, classwise = cw)
  }
  
  train_eval <- function(train_data, test_data, model_id, val_id, held_out) {
    cat(sprintf("   --> %s (%s, held out: %s)\n", model_id, val_id, held_out))
    rec   <- recipes::recipe(target ~ ., data = train_data) %>%
      recipes::step_integer(recipes::all_nominal_predictors())
    safe_v <- max(2, min(5, floor(nrow(train_data) / 10)))
    folds <- rsample::vfold_cv(train_data, v = safe_v, strata = target)
    wkflw <- workflows::workflow() %>% workflows::add_model(xgb_spec) %>% workflows::add_recipe(rec)
    tuned <- tryCatch(
      finetune::tune_race_anova(wkflw, resamples = folds, grid = tuning_grid_size,
                                control = finetune::control_race(save_pred = TRUE, verbose_elim = FALSE),
                                metrics = yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)),
      error = function(e) { message("Tuning failed: ", e$message); NULL })
    if (is.null(tuned)) return(NULL)
    best  <- tune::select_best(tuned, metric = "accuracy")
    fit_m <- parsnip::fit(tune::finalize_workflow(wkflw, best), data = train_data)
    preds <- broom::augment(fit_m, test_data)
    cm    <- table(preds$target, preds$.pred_class)
    lev   <- levels(train_data$target)
    cmc   <- matrix(0, length(lev), length(lev), dimnames = list(lev, lev))
    for (r in rownames(cm)) for (cc in colnames(cm))
      if (r %in% lev && cc %in% lev) cmc[r, cc] <- cm[r, cc]
    extract_matrix_metrics(cmc, model_id, val_id, held_out)
  }
  
  scenarios <- list(With_Species = function(df) df,
                    Without_Species = function(df) df %>% dplyr::select(-species))
  grouping  <- list(
    With_Species    = c("species", "class", "interaction_type", "data_collection", "captive"),
    Without_Species = c("class", "interaction_type", "data_collection", "captive"))
  
  for (scen in names(scenarios)) {
    cat(sprintf("\nScenario: %s\n", scen))
    df_scen <- scenarios[[scen]](df_shuffled)
    
    split <- rsample::initial_split(df_scen, strata = target, prop = 0.7)
    m <- train_eval(rsample::training(split), rsample::testing(split),
                    scen, "Random_70_30_Split", "None")
    if (!is.null(m)) { overall <- dplyr::bind_rows(overall, m$overall)
    classwise <- dplyr::bind_rows(classwise, m$classwise) }
    
    for (col in grouping[[scen]]) {
      counts <- df_shuffled %>% dplyr::count(!!rlang::sym(col)) %>% dplyr::arrange(dplyr::desc(n))
      groups <- as.character(counts[[col]])
      if (col == "species") groups <- utils::head(groups, top_n_species)  # too many levels for exhaustive
      for (grp in groups) {
        tr <- df_scen %>% dplyr::filter(df_shuffled[[col]] != grp)
        te <- df_scen %>% dplyr::filter(df_shuffled[[col]] == grp)
        if (nrow(te) > 0 && length(unique(tr$target)) > 1) {
          mm <- train_eval(tr, te, scen, paste0("Leave_One_", col, "_Out"), grp)
          if (!is.null(mm)) { overall <- dplyr::bind_rows(overall, mm$overall)
          classwise <- dplyr::bind_rows(classwise, mm$classwise) }
        }
      }
    }
  }
  list(overall_performance = overall, classwise_performance = classwise)
}


# =============================================================================
# F. METADATA-MODEL RECOVERY + EVALUATION + LaTeX
# -----------------------------------------------------------------------------
# Run the recover_* steps ONCE in the recovery conda env (older xgboost) to
# export portable boosters; run the evaluate_* steps in the normal env.
# =============================================================================

#' (recovery env) Export a portable booster + trained recipe from a saved fit.
recover_booster <- function(fit_rds_path, booster_out_path, recipe_out_path) {
  wf <- readRDS(fit_rds_path)$.workflow[[1]]
  booster <- xgboost::xgb.Booster.complete(workflows::extract_fit_engine(wf), saveraw = TRUE)
  xgboost::xgb.save(booster, booster_out_path)
  saveRDS(workflows::extract_recipe(wf, estimated = TRUE), recipe_out_path)
  invisible(list(booster = booster, recipe = recipe_out_path))
}

#' Embed feature names + confirm column order into a portable booster.
finalize_booster <- function(booster_json_path, preproc, sample_row, output_path) {
  baked <- recipes::bake(preproc, new_data = sample_row)
  correct_order <- setdiff(names(baked), "target")
  booster <- xgboost::xgb.load(booster_json_path)
  xgboost::setinfo(booster, "feature_name", correct_order)
  xgboost::xgb.save(booster, output_path)
  list(booster = xgboost::xgb.load(output_path), correct_order = correct_order)
}

#' Predict test data and attach the confirmed class labels.
predict_with_labels <- function(new_data, booster, preproc, correct_order, class_order) {
  baked <- recipes::bake(preproc, new_data = new_data)
  probs <- stats::predict(booster, as.matrix(baked %>% dplyr::select(dplyr::all_of(correct_order))))
  colnames(probs) <- class_order
  tibble::tibble(
    truth = new_data$target,
    predicted_class = factor(colnames(probs)[apply(probs, 1, which.max)], levels = class_order),
    tibble::as_tibble(probs))
}

#' Per-class precision/recall/F1 with the NA -> 0 convention (avoids yardstick
#' silently inflating undefined classes).
compute_fixed_metrics <- function(truth, predicted_class) {
  classes <- levels(truth)
  per_class <- purrr::map_df(classes, function(cls) {
    tp <- sum(predicted_class == cls & truth == cls)
    fp <- sum(predicted_class == cls & truth != cls)
    fn <- sum(predicted_class != cls & truth == cls)
    support <- sum(truth == cls)
    prec <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    rec  <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    f1   <- if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
    tibble::tibble(class = cls, precision = prec, recall = rec, f1 = f1, support = support)
  })
  list(per_class = per_class,
       macro    = c(precision = mean(per_class$precision), recall = mean(per_class$recall), f1 = mean(per_class$f1)),
       weighted = c(precision = weighted.mean(per_class$precision, per_class$support),
                    recall    = weighted.mean(per_class$recall,    per_class$support),
                    f1        = weighted.mean(per_class$f1,        per_class$support)))
}

#' Confusion matrix + comprehensive metrics table (model vs majority baseline).
evaluate_model <- function(results) {
  conf_mat <- table(Truth = results$truth, Predicted = results$predicted_class)
  acc <- mean(results$truth == results$predicted_class)
  truth_p <- prop.table(table(results$truth)); pred_p <- prop.table(table(results$predicted_class))
  exp_agree <- sum(truth_p * pred_p[names(truth_p)])
  kappa <- (acc - exp_agree) / (1 - exp_agree)
  
  model_m <- compute_fixed_metrics(results$truth, results$predicted_class)
  majority <- names(which.max(table(results$truth)))
  base_pred <- factor(majority, levels = levels(results$truth))
  base_acc <- mean(results$truth == base_pred)
  base_exp <- sum(truth_p * prop.table(table(base_pred))[names(truth_p)])
  base_kappa <- if (base_exp == 1) 0 else (base_acc - base_exp) / (1 - base_exp)
  base_m <- compute_fixed_metrics(results$truth, base_pred)
  
  metrics_table <- tibble::tibble(
    Metric_Group = c("Discrete Labels", "Discrete Labels", "Macro Averages", "Macro Averages",
                     "Macro Averages", "Weighted Averages", "Weighted Averages", "Weighted Averages"),
    Metric_Name  = c("Overall Accuracy", "Cohen's Kappa", "Macro-Precision", "Macro-Recall",
                     "Macro-F1 Score", "Weighted-Precision", "Weighted-Recall", "Weighted-F1 Score"),
    Model_Score  = c(acc, kappa, model_m$macro["precision"], model_m$macro["recall"], model_m$macro["f1"],
                     model_m$weighted["precision"], model_m$weighted["recall"], model_m$weighted["f1"]),
    Baseline_Score = c(base_acc, base_kappa, base_m$macro["precision"], base_m$macro["recall"], base_m$macro["f1"],
                       base_m$weighted["precision"], base_m$weighted["recall"], base_m$weighted["f1"])
  ) %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4)))
  
  list(conf_mat = conf_mat, metrics_table = metrics_table,
       model_per_class = model_m$per_class, baseline_per_class = base_m$per_class)
}

#' Emit a metrics table as LaTeX (base R, no glue dependency).
make_latex_table <- function(metrics_table, caption, label1, label2, file = NULL) {
  m <- metrics_table$Model_Score; b <- metrics_table$Baseline_Score
  f <- function(x) sprintf("%.4f", x)
  lines <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{\\changeto{}{", caption, "}}"), paste0("\\label{", label1, "}"),
    "\\begin{tabular}{llcc}", "\\toprule",
    "\\textbf{Metric Group} & \\textbf{Metric Name} & \\textbf{Model Score} & \\textbf{Baseline Score} \\\\",
    "\\midrule", "\\textbf{Discrete Labels} ",
    paste0("    & Overall Accuracy   & ", f(m[1]), " & ", f(b[1]), " \\\\"),
    paste0("    & Cohen's Kappa ($\\kappa$) & ", f(m[2]), " & ", f(b[2]), " \\\\"),
    "\\addlinespace", "\\textbf{Macro Averages} \\textit{(Equal Weight)} ",
    paste0("    & Macro-Precision    & ", f(m[3]), " & ", f(b[3]), " \\\\"),
    paste0("    & Macro-Recall       & ", f(m[4]), " & ", f(b[4]), " \\\\"),
    paste0("    & Macro-F1 Score     & ", f(m[5]), " & ", f(b[5]), " \\\\"),
    "\\midrule", "\\textbf{Weighted Averages} \\textit{(Size Proportional)} ",
    paste0("    & Weighted-Precision & ", f(m[6]), " & ", f(b[6]), " \\\\"),
    paste0("    & Weighted-Recall    & ", f(m[7]), " & ", f(b[7]), " \\\\"),
    paste0("    & Weighted-F1 Score  & ", f(m[8]), " & ", f(b[8]), " \\\\"),
    "\\bottomrule", "\\end{tabular}", paste0("\\label{", label2, "}"), "\\end{table}")
  if (!is.null(file)) { writeLines(lines, file); cat("Written to:", file, "\n") } else cat(lines, sep = "\n")
  invisible(paste(lines, collapse = "\n"))
}


# =============================================================================
# G. THRESHOLD-SENSITIVITY TEST (unweighted .edges networks)
# -----------------------------------------------------------------------------
# Does binarising / thresholding change the predicted class?
# The .edges files are unweighted, so we derive a structural connection-strength
# proxy (neighbourhood-overlap / Jaccard, or 1/effective-resistance), delete the
# weakest edges, and re-classify. This is a TOPOLOGICAL robustness check, NOT a
# re-analysis of the original contact structure.
# =============================================================================

read_edges_graph <- function(path, sep = "") {
  raw <- utils::read.table(path, header = FALSE, sep = sep, fill = TRUE,
                           stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(raw) < 2) stop("File ", path, " has fewer than 2 columns.")
  el <- stats::setNames(raw[, 1:2], c("from", "to"))
  el <- el[stats::complete.cases(el), ]
  el$from <- as.character(el$from); el$to <- as.character(el$to)
  el <- el[el$from != el$to, ]
  g <- igraph::graph_from_data_frame(el, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  igraph::delete_vertices(g, which(igraph::degree(g) == 0))
}

largest_component <- function(g) {
  comps <- igraph::components(g, mode = "weak")
  igraph::induced_subgraph(g, igraph::V(g)[comps$membership == which.max(comps$csize)])
}

add_edge_strength <- function(g, weight_method = c("adjacency", "laplacian"),
                              normalize_adjacency = TRUE) {
  weight_method <- match.arg(weight_method)
  A   <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  deg <- rowSums(A)
  el  <- igraph::as_edgelist(g, names = FALSE)
  if (weight_method == "adjacency") {
    CN <- A %*% A
    strength <- vapply(seq_len(nrow(el)), function(k) {
      i <- el[k, 1]; j <- el[k, 2]; inter <- CN[i, j]
      if (normalize_adjacency) { uni <- deg[i] + deg[j] - inter; if (uni > 0) inter / uni else 0 }
      else inter
    }, numeric(1))
  } else {
    Linv <- MASS::ginv(diag(deg) - A)
    strength <- vapply(seq_len(nrow(el)), function(k) {
      i <- el[k, 1]; j <- el[k, 2]
      R <- Linv[i, i] + Linv[j, j] - 2 * Linv[i, j]; if (R > 0) 1 / R else Inf
    }, numeric(1))
  }
  igraph::E(g)$strength <- strength
  g
}

sensitivity_one_graph <- function(g_full, model, graph_name, weight_method, normalize_adjacency,
                                  n_thresholds, use_full_network, use_largest_component,
                                  min_group_size_for_features) {
  g_lcc <- largest_component(g_full)
  if (igraph::vcount(g_lcc) < min_group_size_for_features) return(NULL)
  baseline_class <- classify_graph(g_lcc, model, tag = graph_name)
  if (is.null(baseline_class)) return(NULL)
  
  g_lcc <- add_edge_strength(g_lcc, weight_method, normalize_adjacency)
  w <- igraph::E(g_lcc)$strength; w <- w[is.finite(w)]
  thr <- stats::quantile(w, probs = seq(0, 0.9, length.out = n_thresholds), na.rm = TRUE)
  ql  <- names(thr); keep <- !duplicated(round(thr, 8)); thr <- thr[keep]; ql <- ql[keep]
  
  rows <- list()
  for (m in seq_along(thr)) {
    g_thr <- igraph::delete_edges(g_lcc, which(igraph::E(g_lcc)$strength < thr[m]))
    g_thr <- igraph::delete_vertices(g_thr, which(igraph::degree(g_thr) == 0))
    modes <- c(if (use_full_network) "full_network", if (use_largest_component) "largest_component")
    for (mode in modes) {
      g_eval <- if (mode == "largest_component") largest_component(g_thr) else g_thr
      pred <- if (igraph::vcount(g_eval) < min_group_size_for_features) NA_character_
      else { p <- classify_graph(g_eval, model, tag = paste0(graph_name, "_", ql[m], "_", mode)); if (is.null(p)) NA_character_ else p }
      rows[[length(rows) + 1]] <- tibble::tibble(
        graph_name = graph_name, weight_method = weight_method, component_mode = mode,
        threshold_quant = ql[m], threshold = as.numeric(thr[m]),
        n_nodes = igraph::vcount(g_eval), n_edges = igraph::ecount(g_eval),
        baseline_class = baseline_class, predicted_class = pred, is_baseline = (m == 1),
        matches_baseline = !is.na(pred) & pred == baseline_class)
    }
  }
  dplyr::bind_rows(rows)
}

threshold_sensitivity_test <- function(
    edges_dir = here::here("data", "Networks-Threshold-Test"),
    model = load_main_model(),
    weight_method = c("adjacency", "laplacian"), normalize_adjacency = TRUE,
    n_thresholds = 10, use_full_network = TRUE, use_largest_component = TRUE,
    min_group_size_for_features = 5, sep = "",
    out_prefix = here::here("outputs", "threshold_sensitivity")) {
  
  weight_method <- match.arg(weight_method)
  files <- list.files(edges_dir, pattern = "\\.edges$", full.names = TRUE)
  if (length(files) == 0) stop("No .edges files in '", edges_dir, "'.")
  
  all_rows <- list()
  for (f in files) {
    gname <- tools::file_path_sans_ext(basename(f)); message("Processing: ", gname)
    g <- tryCatch(read_edges_graph(f, sep = sep), error = function(e) NULL)
    if (is.null(g)) next
    res <- tryCatch(sensitivity_one_graph(g, model, gname, weight_method, normalize_adjacency,
                                          n_thresholds, use_full_network, use_largest_component, min_group_size_for_features),
                    error = function(e) { message("  error: ", conditionMessage(e)); NULL })
    if (!is.null(res)) all_rows[[gname]] <- res
  }
  results <- dplyr::bind_rows(all_rows)
  if (nrow(results) == 0) { warning("No classifiable networks."); return(invisible(list())) }
  
  by_graph <- results %>%
    dplyr::group_by(graph_name, weight_method, component_mode, baseline_class) %>%
    dplyr::summarise(
      n_thresholds_evaluated = dplyr::n(),
      n_distinct_predictions = dplyr::n_distinct(predicted_class[!is.na(predicted_class)]),
      classes_seen = paste(sort(unique(stats::na.omit(predicted_class))), collapse = ", "),
      prop_matching_baseline = mean(matches_baseline),
      stable = all(matches_baseline[!is.na(predicted_class)]),
      n_unclassifiable = sum(is.na(predicted_class)), .groups = "drop") %>%
    dplyr::arrange(component_mode, graph_name)
  wide <- results %>%
    dplyr::select(graph_name, component_mode, threshold_quant, predicted_class) %>%
    tidyr::pivot_wider(names_from = threshold_quant, values_from = predicted_class)
  
  utils::write.csv(results,  paste0(out_prefix, "_raw.csv"),      row.names = FALSE)
  utils::write.csv(by_graph, paste0(out_prefix, "_by_graph.csv"), row.names = FALSE)
  utils::write.csv(wide,     paste0(out_prefix, "_wide.csv"),     row.names = FALSE)
  invisible(list(results = results, by_graph = by_graph, wide = wide))
}


# =============================================================================
# H. PLOT HELPERS
# =============================================================================

#' Base theme shared by figures.
theme_paper <- function(base_size = 18) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "serif", size = base_size),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom")
}

#' Stacked proportion bars of predicted class, faceted by taxonomic class.
metadata_plot_taxa <- function(combined_df, facet_column = "Class",
                               group_by_column = "Data_collection") {
  d <- combined_df %>%
    dplyr::group_by(Class, !!rlang::sym(group_by_column), Predicted_classes) %>%
    dplyr::mutate(Count = dplyr::n()) %>%
    dplyr::group_by(Class, !!rlang::sym(group_by_column)) %>%
    dplyr::mutate(Proportion = Count / sum(Count)) %>% dplyr::ungroup()
  x_name <- tools::toTitleCase(gsub("_", " ", group_by_column))
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[group_by_column]], y = Proportion,
                                  fill = Predicted_classes)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::facet_grid(stats::as.formula(paste("~", facet_column)), scales = "free_x") +
    ggplot2::labs(x = x_name, y = "Proportion", fill = "Predicted Classes") +
    ggplot2::scale_fill_manual(values = COLOR_PALETTE) +
    theme_paper() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1))
}

#' Multiclass AUC for a predict_with_labels() result.
#' `results` must have: truth, predicted_class, and one probability column per class.
compute_auc <- function(results, class_order = CLASS_ORDER) {
  
  res <- results %>%
    dplyr::mutate(truth = factor(as.character(truth), levels = class_order))
  
  stopifnot(all(class_order %in% names(res)))
  stopifnot(all(abs(rowSums(res[, class_order]) - 1) < 1e-6))
  
  ## ---- Overall (multiclass) --------------------------------------------------
  overall <- dplyr::bind_rows(
    yardstick::roc_auc(res, truth = truth, dplyr::all_of(class_order), estimator = "hand_till"),
    yardstick::roc_auc(res, truth = truth, dplyr::all_of(class_order), estimator = "macro"),
    yardstick::roc_auc(res, truth = truth, dplyr::all_of(class_order), estimator = "macro_weighted")
  ) %>%
    dplyr::select(estimator = .estimator, auc = .estimate)
  
  ## ---- Class-wise one-vs-rest ------------------------------------------------
  per_class <- purrr::map_dfr(class_order, function(cls) {
    truth_bin <- factor(ifelse(as.character(res$truth) == cls, cls, "other"),
                        levels = c(cls, "other"))
    prob_col  <- res[[cls]]
    support   <- sum(as.character(res$truth) == cls)
    
    if (support == 0 || support == nrow(res)) {   # AUC undefined
      return(tibble::tibble(class = cls, support = support,
                            roc_auc = NA_real_, pr_auc = NA_real_))
    }
    tibble::tibble(
      class   = cls,
      support = support,
      roc_auc = yardstick::roc_auc_vec(truth_bin, prob_col, event_level = "first"),
      pr_auc  = yardstick::pr_auc_vec(truth_bin,  prob_col, event_level = "first")
    )
  })
  
  list(overall = overall, per_class = per_class)
}

#' Overall accuracy from a predict_with_labels() result.
overall_accuracy <- function(results) {
  mean(as.character(results$truth) == as.character(results$predicted_class))
}

#' SHAP variable-importance bar for one class.
shap_importance_plot <- function(shap_class) {
  shapviz::sv_importance(shap_class, kind = "bar", bee_width = 0.3) +
    ggplot2::ggtitle(" ") + ggplot2::xlab("Mean SHAP Value") + theme_paper(24)
}