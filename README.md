# Animal Social Networks and Standard Theoretical Generative Models

Code and data for the manuscript *"Animal Social Networks Most Closely Resemble
Geometric Graph Models Among Standard Generative Models"* (Appaw, Silk, Rushmore,
VanderWaal, Charleston & Fountain-Jones).

We classify empirical animal social networks from the Animal Social Network
Repository (ASNR) against five standard theoretical generative models:
Erdős–Rényi (`ER`), stochastic block model (`sbm`), scale-free (`SF`),
geometric/spatial (`Spatial`) and small-world (`SW`). We use a two-stage
interpretable XGBoost pipeline, then ask which metadata (species, interaction
type, data-collection method, captivity status, observation duration) are
associated with the structure a network is assigned.

**Note on interpretation.** A predicted class means *structural resemblance to
that model within the restricted, pre-specified set of five candidates*. It is
not evidence that the corresponding process generated the network. SHAP values
are likewise attributions within the fitted models, not causal effects on
network structure.

---

## Quick start

```bash
git clone https://github.com/araimacarol/Animal-Social-Networks-Theoretical-Generative-Models.git
cd Animal-Social-Networks-Theoretical-Generative-Models
git lfs pull          # downloads models/final.fitted.model.rds (see "Large file")
```

```r
install.packages("renv"); renv::restore()   # exact package versions
source("R/02_first_stage_and_metadata.R")   # quick smoke test
source("run_all.R")                         # whole pipeline (or: Rscript run_all.R 05 06)
```

All paths are project-relative via `{here}`. There is no `setwd()` in the
pipeline, and the working directory you launch from does not matter.

---

## Repository layout

```
.
├── run_all.R                       # runs the pipeline (steps 02-09, 11, PCA)
├── R/
│   ├── 00_setup.R                  # packages, paths, global constants
│   ├── 01_functions.R              # all reusable functions (no side effects)
│   ├── 02_first_stage_and_metadata.R
│   ├── 03_metadata_model_evaluation.R
│   ├── 04_grouped_cross_validation.R
│   ├── 05_group_level_community_detection.R
│   ├── 06_threshold_sensitivity.R
│   ├── 07_interpretation_shap_hstats.R
│   ├── 08_sbm_simulation_tables.R
│   ├── 09_group_level_tables.R
│   ├── 10_first_stage_test_performance.R
│   ├── 11_illustrative_figures.R
│   ├── PCA_Applicability.R
│   ├── point2_export_derived.R     # one-off: large training files -> small derived files
│   └── sbm-sim-different-communities.R
├── data/                           # analysis inputs
│   ├── Empirical-Network-Edges/    # empirical .edges networks
│   └── Networks-Threshold-Test/    # networks for the edge-removal test (step 06)
├── models/                         # recovered models, data splits, derived files
└── outputs/                        # generated CSV tables (and inputs to 06 and 08)
```

`figures/` and `tables/` are created at run time and are not tracked.

---

## Scripts, and the tables and figures each produces

Table and figure numbers refer to the revised manuscript (main text) and its
electronic supplementary material (S).

| Step | Script | Produces | Runtime |
|------|--------|----------|---------|
| 00 | `00_setup.R` | packages, paths, constants (sourced by every script) | instant |
| 01 | `01_functions.R` | all functions (sourced, never run alone) | instant |
| 02 | `02_first_stage_and_metadata.R` | first-stage class of every empirical network; Figure 3, Fig. S2; Tables S2–S5 | minutes |
| 03 | `03_metadata_model_evaluation.R` | Tables 5, 6, S7 | minutes |
| 04 | `04_grouped_cross_validation.R` | Tables 7, S8 (both models refitted on identical splits) | **slow** |
| 05 | `05_group_level_community_detection.R` | Table S6; cross-method results in Section 3.3 | moderate |
| 06 | `06_threshold_sensitivity.R` | Tables 3, S1, S23, S24 | slow if `RERUN = TRUE` |
| 07 | `07_interpretation_shap_hstats.R` | Figures 4, 5; Figs S3–S11 | moderate |
| 08 | `08_sbm_simulation_tables.R` | Tables S16, S17 | fast |
| 09 | `09_group_level_tables.R` | Table 4; Tables S12–S14 | moderate |
| 10 | `10_first_stage_test_performance.R` | Table 2 | moderate |
| 11 | `11_illustrative_figures.R` | Figure 1, Fig. S12 (illustrations, one random draw each) | fast |
| PCA | `PCA_Applicability.R` | Figure 2, Fig. S1; Tables S18–S22, S25 | moderate |
| — | `sbm-sim-different-communities.R` | the SBM simulations behind Tables S16/S17 | **slow** |
| — | `point2_export_derived.R` | the three derived files read by `PCA_Applicability.R` | moderate |

`run_all.R` runs 02–09, 11 and PCA. Run these separately:
- `sbm-sim-different-communities.R`. A rerun gives new random networks from the
  same design, because the original run set no seed. It writes to `*_rerun.csv`
  so that the deposited simulation output is not overwritten.
- `point2_export_derived.R`. See "Large file" below.
- `10_first_stage_test_performance.R`. It needs the held-out simulated test set,
  `models/df_test_new.rds`.

Table 1 (simulation parameters) and Tables S9–S11 and S15 (definitions and the
species list) are static and are not produced by code.

---

## Data files (`data/`)

| File | Description |
|------|-------------|
| `GraphFeatOnAllAnimNets.csv` | Graph features for every empirical ASNR network |
| `Network_repository_metaData_combined.csv` | Network metadata (species, class, interaction type, data collection, captivity, duration) |
| `Badger_Groups.csv` | Field-defined badger social groups (`NodeID`, `Group`) |
| `Empirical-Network-Edges/*.edges` | Empirical edge lists (unweighted) |
| `Networks-Threshold-Test/*.edges` | Networks for the edge-removal test |

Some edge files keep a `.csv .edges` suffix from their original download.
`resolve_edge_file()` in `01_functions.R` handles this, so nothing needs to be
renamed.

## Model and derived files (`models/`)

| File | Used by | Description |
|------|---------|-------------|
| `booster_portable_mainmodel.json` + `preproc_trained_mainmodel.rds` | 02, 05, 06, 08, 09, 10, PCA, SBM simulations | **Recovered** first-stage classifier: every first-stage classification |
| `booster_portable.json` + `preproc_trained_meta1.rds` | 03, 07 | Recovered metadata Model 1 (with species) |
| `booster_portable_meta2.json` + `preproc_trained_meta2.rds` | 03, 07 | Recovered metadata Model 2 (without species) |
| `metadata_train{,2}.rds`, `metadata_test{,2}.rds` | 03, 07 | Train/test splits of the original metadata fits |
| `df_final_meta_data.rds` | 04, 07 | 636 networks: metadata and first-stage class |
| `final.emp.data.rds` | PCA | 636 empirical networks, Figure 2 feature set |
| `point2_training_size_by_class.csv`, `point2_training_feature_ranges.csv`, `point2_fig2_pca.rds` | PCA | Derived from the training data by `point2_export_derived.R` |
| `df_train_new.rds` | `point2_export_derived.R` | 309,499 pre-SMOTE training networks |
| `final.fitted.model.rds` (Git LFS) | `point2_export_derived.R` | Original fitted workflow. Only its stored training predictors are read, for the Figure 2 PCA; its booster is never used to predict |
| `df_test_new.rds` | 10 | Held-out simulated test set (132,643 networks) |

### Inputs stored in `outputs/`

Two scripts summarise stored results rather than recomputing them:
- `06_threshold_sensitivity.R` reads `threshold_sensitivity_raw.csv`. Set
  `RERUN <- TRUE` to recompute it.
- `08_sbm_simulation_tables.R` reads
  `all_sbm_preds_{lower,higher}_between_connectn_probs.csv`. It first checks every
  row (mean degree = 2E/N) and that the recovered model reproduces every stored
  class.

---

## Large file: `final.fitted.model.rds`

This 206 MB file is stored with **Git LFS** (`git lfs pull` after cloning). It is
needed only to rebuild the three derived files with `point2_export_derived.R`.
`PCA_Applicability.R` uses the derived files when they are present, and
otherwise reads `df_train_new.rds` and `final.fitted.model.rds` directly. Both
routes give the same results: PC1 38.3%, PC2 22.3%, 610 of 636 empirical
networks inside the convex hull (95.9%).

---

## Model recovery

The original tidymodels workflow objects fail to reload under current XGBoost
(`object is corrupted or is from an incompatible XGBoost version`). The
repository therefore ships **version-portable JSON boosters** plus the **trained
recipes**, and predicts with the native booster (`xgboost::xgb.load()` +
`recipes::bake()`). Every classification in the paper uses these recovered
models, through `load_main_model()` / `predict_main_model()` for the first-stage
classifier.

To re-export the boosters, run `recover_booster()` once in an environment with
the older XGBoost that wrote the original files. The calls are commented out at
the top of `03_metadata_model_evaluation.R`.

---

## Notes on specific analyses

**Group-level classification (05, 09).** Communities smaller than **8 nodes**
are not classified (`MIN_GROUP_SIZE` in `00_setup.R`). Badger groups are the
field-defined groups; chimpanzee, giraffe and ant groups are detected with
walktrap. Network modularity is that of the walktrap partition of the whole
network. The cross-method comparison (05) uses leading eigenvector, Louvain and
walktrap. Infomap returned a single community for every chimpanzee, giraffe and
ant network, so its rows are dropped from all outputs. It is still run, because
Louvain shares its random-number stream: removing Infomap would change the
Louvain partitions.

**Edge-removal test (06).** The networks are unweighted, so observed interaction
weights cannot be thresholded. Edges are ranked by neighbourhood overlap
(Jaccard index of shared neighbours) and the weakest are removed in steps of 10
percentiles. This is a **topological perturbation test**. It does not test
sensitivity to interaction weights and does not show that binarisation is
inconsequential.

**Metadata models (03, 04).** Tables 5 and 6 evaluate the original fits, each on
its own test set, so they are not used to compare the two models. Script 04
refits both models on identical splits: one stratified 70:30 random split and 27
grouped hold-outs. Each split reports a majority-class baseline taken from its
training fold.

**SHAP and H-statistics (07).** SHAP values are computed for each network in the
training sets of the recovered metadata models (Tables 5, 6, S7) and are on the
log-odds (margin) scale. Dependence plots show one point per network, with the
mean of each level as a black diamond; nothing is summed. H-statistics are
unnormalised, on the predicted-probability scale.

---

## Environment

- R ≥ 4.2.3
- Exact package versions are pinned in `renv.lock`; `renv::restore()` reproduces
  them.
- `_dependencies.R` is a manifest read only by `renv` and is never sourced.
  `00_setup.R` loads packages with `pacman::p_load()`, which `renv` cannot parse.
  If you add a package to `p_load()`, add it to `_dependencies.R` too.

## Citation

Appaw, R.C., Silk, M.J., Rushmore, J., VanderWaal, K., Charleston, M.A. &
Fountain-Jones, N.M. *Animal Social Networks Most Closely Resemble Geometric
Graph Models Among Standard Generative Models.*

The first-stage classifier is adapted from: Appaw, R.C., Fountain-Jones, N.M. &
Charleston, M.A. (2025) *Leveraging advances in machine learning for the robust
classification and interpretation of networks.* Royal Society Open Science
12:240458. <https://doi.org/10.1098/rsos.240458>

## License

MIT (code). ASNR network data are redistributed under their original terms; see
Sah et al. (2019), *Scientific Data* 6:44, and the Animal Social Network
Repository.
