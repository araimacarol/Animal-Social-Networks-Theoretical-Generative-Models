# =============================================================================
# _dependencies.R  --  package manifest for {renv}
# -----------------------------------------------------------------------------
# This file is NEVER sourced by the pipeline. It exists solely so that
# renv::snapshot() can discover every package the project needs.
#
# renv finds dependencies by scanning for library(), require() and pkg:: calls.
# It cannot see inside pacman::p_load(), which is how 00_setup.R loads packages,
# so without this file renv.lock would silently omit ~13 packages (tidyverse,
# themis, hardhat, missForest and others).
#
# If you add a package to p_load() in 00_setup.R, add it here too.
# =============================================================================

library(Boruta)
library(FactoMineR)
library(MASS)
library(Matrix)
library(RColorBrewer)
library(broom)
library(cowplot)
library(dplyr)
library(factoextra)
library(finetune)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(hardhat)
library(here)
library(hstats)
library(igraph)
library(janitor)
library(jsonlite)
library(knitr)
library(missForest)
library(pacman)
library(parsnip)
library(patchwork)
library(purrr)
library(recipes)
library(renv)
library(rlang)
library(rsample)
library(shapviz)
library(themis)
library(tibble)
library(tidymodels)
library(tidyr)
library(tidyverse)
library(tune)
library(workflows)
library(xgboost)
library(yardstick)
