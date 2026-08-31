# =============================================================================
# 00_setup.R  --  packages, paths, and global options
# -----------------------------------------------------------------------------
# Source this first (run_all.R does this for you). It loads every package the
# pipeline needs and defines project-relative paths via {here}, so no machine-
# specific setwd() is required.
# =============================================================================

## ---- Packages ---------------------------------------------------------------
## Installed on demand with {pacman}; versions are pinned in the README /
## renv.lock. Grouped by role for readability.
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  # core / data wrangling
  tidyverse, janitor, here, jsonlite,
  # graphs & features
  igraph, Matrix,
  # modelling (tidymodels stack + xgboost)
  tidymodels, recipes, workflows, parsnip, tune, finetune, themis, hardhat,
  xgboost,
  # feature selection & imputation
  Boruta, missForest,
  # interpretation
  shapviz, hstats,
  # dimension reduction / PCA
  FactoMineR, factoextra,
  # plotting helpers
  ggpubr, patchwork, cowplot, gridExtra, RColorBrewer,
  # misc
  knitr, MASS
)

## {MASS} masks dplyr::select -- keep dplyr's version project-wide.
select <- dplyr::select

## ---- Paths ------------------------------------------------------------------
## Everything is relative to the project root (the folder containing the .Rproj
## / .git). Create the standard sub-folders if they do not yet exist.
dir_data     <- here::here("data")        # raw inputs (csv, .edges, Badger_Groups.csv)
dir_models   <- here::here("models")      # trained/recovered model objects (.json, .rds)
dir_outputs  <- here::here("outputs")     # generated csv tables

## Empirical edge lists. Keep the ASNR/author file names as downloaded; the
## helper resolve_edge_file() (01_functions.R) copes with the ".csv .edges"
## suffix some of them carry, so nothing has to be renamed by hand.
dir_edges    <- here::here("data", "Empirical-Network-Edges")
dir_thresh   <- here::here("data", "Networks-Threshold-Test")

## Scratch folders for generated artefacts. These are created at run time and
## are NOT tracked in git (see .gitignore): the repository ships the CSVs in
## outputs/, while the manuscript builds its own figures and formatted tables.
dir_figures  <- here::here("figures")     # generated figures (untracked)
dir_tables   <- here::here("tables")      # generated LaTeX tables (untracked)

for (d in c(dir_data, dir_models, dir_outputs, dir_figures, dir_tables)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## ---- Global constants -------------------------------------------------------
## Class labels used throughout. The pipeline stores the compact codes; display
## code maps Spatial -> "Geometric Graphs".
CLASS_ORDER   <- c("ER", "sbm", "SF", "Spatial", "SW")
CLASS_DISPLAY <- c(
  ER = "Erdos-Renyi", sbm = "Stochastic-Block-Model",
  SF = "Scale-Free", Spatial = "Geometric Graphs", SW = "Small-World"
)

## Okabe-Ito-based palette (ER, Geometric, Scale-free, Small-world, SBM)
COLOR_PALETTE <- c("#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2")

## Minimum community size for group-level classification (Reviewer point:
## do not classify tiny groups). Field-defined badger groups span 2-9 nodes.
MIN_GROUP_SIZE <- 8

set.seed(123)
