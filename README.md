# Animal Social Networks and Standard Theoretical Generative Models

Code and data for the manuscript *"The commonality of Geometric Graphs in Animal
Social Networks"* (Appaw, Silk, Rushmore,VanderWaal, Charleston & Fountain-Jones).

We classify empirical animal social networks from the Animal Social Network
Repository (ASNR) against five standard theoretical generative models —
Erdős–Rényi (`ER`), stochastic block model (`sbm`), scale-free (`SF`),
geometric/spatial (`Spatial`) and small-world (`SW`) — using a two-stage
interpretable XGBoost pipeline, then ask which metadata (species, interaction
type, data-collection method, captivity status, observation duration) are
associated with the structure a network is assigned.

**Note on interpretation.** Throughout the code and the paper, a predicted class
means *structural resemblance to that model within the restricted, pre-specified
set of five candidates* — not evidence that the corresponding process generated
the network. Shapley Additive Explanations (SHAP) values are likewise attributions within the fitted models, not
causal effects on network structure.
---

## Quick start

```r
# Clone and open the project (RStudio: open the .Rproj — no setwd() needed)
git clone https://github.com/araimacarol/Animal-Social-Networks-Theoretical-Generative-Models.git
cd Animal-Social-Networks-Theoretical-Generative-Models

#Restore the exact package versions
install.packages("renv"); renv::restore()

#Run the whole pipeline
Rscript run_all.R
```

`setup_migrate_files.R` is a **one-time helper** used when this repository was
first assembled from the original flat working folder. 
---

## Repository layout

```
.
├── run_all.R                            # master script: runs the pipeline end to end
├── setup_migrate_files.R                # one-time migration helper (not needed for a clone)
├── R/
│   ├── 00_setup.R                       # packages, paths, global constants
│   ├── 01_functions.R                   # ALL reusable functions (no side effects)
│   ├── 02_first_stage_and_metadata.R
│   ├── 03_metadata_model_evaluation.R
│   ├── 04_grouped_cross_validation.R
│   ├── 05_group_level_community_detection.R
│   ├── 06_threshold_sensitivity.R
│   ├── 07_interpretation_shap_hstats.R
│   └── sbm-sim-different-communities.R  # standalone SBM simulation study (not in run_all.R)
├── data/                                # analysis inputs
│   ├── Empirical-Network-Edges/         # 848 empirical .edges networks
│   └── Networks-Threshold-Test/         # 350 networks for the threshold check
├── models/                              # trained + recovered model objects
└── outputs/                             # generated CSV tables
```

## Script order

| Step | Script | What it does | Runtime |
|------|--------|--------------|---------|
| 00 | `00_setup.R` | Loads packages, defines paths/constants | instant |
| 01 | `01_functions.R` | Defines all functions (sourced, never run alone) | instant |
| 02 | `02_first_stage_and_metadata.R` | Classifies empirical networks; joins metadata; writes per-category sample sizes | minutes |
| 03 | `03_metadata_model_evaluation.R` | Both metadata models vs majority-class baseline (macro-F1, balanced accuracy, class-wise P/R) | minutes |
| 04 | `04_grouped_cross_validation.R` | Leave-one-group-out CV (species/class/interaction/collection/captivity) | **slow** |
| 05 | `05_group_level_community_detection.R` | Group-level classification, min group size 8, across 4 detection algorithms | moderate |
| 06 | `06_threshold_sensitivity.R` | Predicted-class robustness to pruning weak ties | **slow** |
| 07 | `07_interpretation_shap_hstats.R` | SHAP + H-statistics | moderate |
| — | `sbm-sim-different-communities.R` | SBM simulations across k = 2–6 and low/high between-community connectivity (Supp. Tables S16/S17). Run separately, after `00_setup.R`. | **slow** |

---

## Data files

These are analysis inputs, not initially generated:

| File | Description |
|------|-------------|
| `data/GraphFeatOnAllAnimNets.csv` | Graph features computed for every empirical ASNR network |
| `data/Network_repository_metaData_combined.csv` | Network metadata (species, class, interaction type, data collection, captivity, duration) |
| `data/Badger_Groups.csv` | Field-defined badger social group memberships (`NodeID`, `Group`) |
| `data/Empirical-Network-Edges/*.edges` | Empirical edge lists — unweighted, whitespace-separated |
| `data/Networks-Threshold-Test/*.edges` | Networks used for the threshold-sensitivity check |

### A note on edge-file names

The edge lists are kept under their **original download names**. Some carry a
`.csv .edges` suffix (with a space) rather than plain `.edges`. Nothing needs to
be renamed: `resolve_edge_file()` in `01_functions.R` tries `<name>.edges`,
`<name>.csv .edges` and `<name>.csv.edges` in turn, so the pipeline finds the
file either way. If a network genuinely cannot be found, step 05 stops with an
explicit message naming it rather than silently dropping it.

### Generated outputs (`outputs/`)

| File | Produced by | Contents |
|------|-------------|----------|
| `Predicted-Class-of-Animal-Social-network.csv` | 02 | Predicted generative model per empirical network |
| `meta_combined_df.csv` | 02 | Predictions joined to metadata |
| `sample_sizes_by_{class,interaction,collection,captive}.csv` | 02 | Networks per metadata category (Supp. Tables S2–S5) |
| `metrics_meta1.csv`, `metrics_meta2.csv` | 03 | Model vs baseline metrics (Tables 6 & 7) |
| `classwise_meta1.csv`, `classwise_meta2.csv` | 03 | Class-wise precision/recall/F1 (Supp. Table S7) |
| `pipeline_overall_performance2.csv` | 04 | Random vs grouped-CV performance (Table 8, Supp. S8) |
| `pipeline_classwise_performance2.csv` | 04 | Class-wise metrics per held-out group |
| `community_summary.csv` | 05 | Group sizes detected per network per method |
| `community_predictions.csv` | 05 | Predicted class per community |
| `merged_community_predictions.csv` | 05 | Merged group-level table (Table 5, Supp. S6, S12–S14) |
| `community_method_stability.csv` | 05 | % geometric per detection algorithm |
| `threshold_sensitivity_{raw,by_graph,wide}.csv` | 06 | Class stability under edge pruning |
| `threshold_sensitivity_by_class.csv` | 06 | Headline invariance by class (Table 3) |
| `metadata_sensitivity_{overall,classwise,auc}.csv` | supporting | Metadata-model sensitivity checks |

### Model objects (`models/`)

| File | Description |
|------|-------------|
| `booster_portable_mainmodel.json` | **Recovered** first-stage classifier (version-portable) |
| `preproc_trained_mainmodel.rds` | Trained recipe for the first-stage classifier |
| `booster_portable.json` / `booster_portable_meta2.json` | Recovered metadata models (with / without species) |
| `booster_portable_named.json` / `booster_portable_named_meta2.json` | As above, with feature names retained |
| `preproc_trained_meta1.rds` / `preproc_trained_meta2.rds` | Their trained recipes |
| `df_final_meta_data.rds` | Imputed metadata + target used for modelling |
| `metadata_train{,2}.rds`, `metadata_test{,2}.rds` | Train/test splits (with / without species) |
| `metadata.prep.new{,2}.rds` | Prepped recipes |
| `xgb.final.metadata.fit.new{,2}.rds` | Original tidymodels metadata workflow fits (see below) |
| `final_meta_shap_{species,no_species}.rds` | Cached SHAP objects (step 07) |
| `meta_interactn_hstats_no_species.rds` | Cached H-statistics (step 07) |

---

## Model recovery (important)

The metadata models were originally saved as tidymodels workflow objects. After
an XGBoost version bump these fail to reload:

```
Error: object is corrupted or is from an incompatible XGBoost version
```

We therefore ship **version-portable JSON boosters** plus the **trained recipes**,
and predict with the native booster (`xgboost::xgb.load()` + `recipes::bake()`).
This is what `load_main_model()` and `03_metadata_model_evaluation.R` use, and it
is why the repository does not depend on a specific XGBoost build.

To re-export the portable boosters, run `recover_booster()` **once** in
an environment with the *older* XGBoost that wrote the original files (the calls
are at the top of `03_metadata_model_evaluation.R`, commented out).

Every script that needs the first-stage classifier — 02, 05, 06 and
`sbm-sim-different-communities.R` — goes through `load_main_model()` /
`predict_main_model()`. 
---

## Notes on two analyses

**Minimum group size (step 05).** Communities smaller than **8 nodes** are not
classified. Field-defined badger social groups span 2–9 individuals, so a
threshold of 8 keeps the largest genuine groups and excludes communities too
small to yield stable structural features. Networks with no qualifying community
are reported as such rather than being force-classified.

**Threshold sensitivity (step 06).** The empirical network `.edges` files are
**unweighted**, so we cannot threshold true interaction weights. We instead remove
edges by a structural connection-strength proxy (neighbourhood overlap / Jaccard
of shared neighbours) and re-classify. This is a **topological robustness check**
— it is not a re-analysis of contact frequency or duration, and should not be
described as one.

---

## Environment

- R ≥ 4.2.3
- Key packages: `igraph`, `tidymodels` (`parsnip`, `recipes`, `workflows`, `tune`,
  `finetune`, `rsample`, `yardstick`), `xgboost`, `shapviz`, `hstats`, `Boruta`,
  `missForest`, `janitor`, `here`, `tidyverse`

Exact versions are pinned in `renv.lock`; `renv::restore()` reproduces them.

`_dependencies.R` is a manifest read only by `renv`. It is never sourced by the
pipeline. It exists because `00_setup.R` loads packages through
`pacman::p_load()`, which `renv` cannot parse — without it, `renv::snapshot()`
would omit around a third of the dependencies. If you add a package to
`p_load()`, add it to `_dependencies.R` as well.

## Citation

Appaw, R.C., Silk, M.J., Rushmore, J., VanderWaal, K., Charleston, M.A. &
Fountain-Jones, N.M. *The commonality of Geometric Graphs in Animal 
Social Networks.*

The first-stage classifier is adapted from: Appaw, R.C., Fountain-Jones, N.M. &
Charleston, M.A. (2025) *Leveraging advances in machine learning for the robust
classification and interpretation of networks.* Royal Society Open Science 12:240458.
<https://doi.org/10.1098/rsos.240458>

## License

MIT (code). ASNR network data are redistributed under their original terms —
see Sah et al. (2019), *Scientific Data* 6:44, and the Animal Social Network
Repository.
