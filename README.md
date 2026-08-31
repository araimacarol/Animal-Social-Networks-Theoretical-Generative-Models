# Animal Social Networks and Standard Theoretical Generative Models

Code and data for the manuscript *"The commonality of Geometric Graphs in Animal
Social Networks"* (Appaw, Silk, Rushmore, VanderWaal, Charleston & Fountain-Jones).

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
the network. Shapley Additive Explanations (SHAP) values are likewise
attributions within the fitted models, not causal effects on network structure.

---

## Quick start

```r
# Clone and open the project (RStudio: open the .Rproj — no setwd() needed)
git clone https://github.com/araimacarol/Animal-Social-Networks-Theoretical-Generative-Models.git
cd Animal-Social-Networks-Theoretical-Generative-Models

# Restore the exact package versions
install.packages("renv"); renv::restore()

# Run the whole pipeline
Rscript run_all.R
```

`setup_migrate_files.R` is a **one-time helper** used when this repository was
first assembled from the original flat working folder.

### If paths do not resolve

Every path is built with `{here}`, which anchors to the folder containing the
`.Rproj` / `.git` / `.here` marker. `here` caches that root the first time it is
called in a session, so if R was started outside the repository the whole
pipeline will look in the wrong place. Check before running anything:

```r
here::here()   # must print the repository root
```

If it does not, open the `.Rproj`, or run `here::i_am("R/00_setup.R")` once from
inside the project and restart R.

---

## Repository layout

```
.
├── run_all.R                            # master script: runs the pipeline end to end
├── setup_migrate_files.R                # one-time migration helper (not needed for a clone)
├── _dependencies.R                      # package manifest for renv (never sourced)
├── renv.lock                            # pinned package versions
├── LICENSE                              # MIT (code)
├── R/
│   ├── 00_setup.R                       # packages, paths, global constants
│   ├── 01_functions.R                   # ALL reusable functions (no side effects)
│   ├── 02_first_stage_and_metadata.R
│   ├── 03_metadata_model_evaluation.R
│   ├── 04_grouped_cross_validation.R
│   ├── 05_group_level_community_detection.R
│   ├── 06_threshold_sensitivity.R
│   ├── 07_interpretation_shap_hstats.R
│   ├── 08_ant_table_s14.R
│   └── sbm-sim-different-communities.R  # standalone SBM simulation study (not in run_all.R)
├── data/                                # analysis inputs
│   ├── Empirical-Network-Edges/         # 852 empirical .edges networks
│   └── Networks-Threshold-Test/         # networks for the threshold check
├── models/                              # trained + recovered model objects
└── outputs/                             # generated CSV tables
```

`figures/` and `tables/` are generated at runtime and are not tracked; the
manuscript builds its own. `outputs/*.csv` **is** tracked.

## Script order

| Step | Script | What it does | Runtime |
|------|--------|--------------|---------|
| 00 | `00_setup.R` | Loads packages, defines paths/constants | instant |
| 01 | `01_functions.R` | Defines all functions (sourced, never run alone) | instant |
| 02 | `02_first_stage_and_metadata.R` | Classifies empirical networks; joins metadata; writes per-category sample sizes | minutes |
| 03 | `03_metadata_model_evaluation.R` | Both metadata models vs majority-class baseline (macro-F1, balanced accuracy, class-wise P/R) | minutes |
| 04 | `04_grouped_cross_validation.R` | Leave-one-group-out CV (species/class/interaction/collection/captivity) | **slow** |
| 05 | `05_group_level_community_detection.R` | Group-level classification, min group size 8, across 4 detection algorithms (chimpanzee, giraffe, ant). **Badgers are not part of the algorithm comparison** — their field-defined groups are classified separately | moderate |
| 06 | `06_threshold_sensitivity.R` | Predicted-class robustness to pruning weak ties | **slow** |
| 07 | `07_interpretation_shap_hstats.R` | SHAP + H-statistics | moderate |
| 08 | `08_ant_table_s14.R` | Rebuilds the ant network/group table (Supp. Table S14) from the recovered model | minutes |
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
file either way. If a network genuinely cannot be found, `resolve_edge_files()`
stops with an explicit message naming it rather than silently dropping it.

### Generated outputs (`outputs/`)

| File | Produced by | Contents |
|------|-------------|----------|
| `Predicted-Class-of-Animal-Social-network.csv` | 02 | Predicted generative model per empirical network |
| `meta_combined_df.csv` | 02 | Predictions joined to metadata |
| `sample_sizes_by_{class,interaction,collection,captive}.csv` | 02 | Networks per metadata category (Supp. Tables S2–S5) |
| `metrics_meta1.csv`, `metrics_meta2.csv` | 03 | Model vs baseline metrics (Tables 5 & 6) |
| `classwise_meta1.csv`, `classwise_meta2.csv` | 03 | Class-wise precision/recall/F1 (Supp. Table S7) |
| `pipeline_overall_performance2.csv` | 04 | Random vs grouped-CV performance (Table 7, Supp. Table S8) |
| `pipeline_classwise_performance2.csv` | 04 | Class-wise metrics per held-out group |
| `community_summary.csv` | 05 | Group sizes detected per network per method |
| `community_predictions.csv` | 05 | Predicted class per detected group |
| `merged_community_predictions.csv` | 05 | Detected groups (chimpanzee, giraffe, ant) plus badger prescribed groups, in one table (Supp. Tables S6 and S13) |
| `community_method_stability.csv` | 05 | % geometric per detection algorithm, for all networks and for networks with modularity ≥ 0.01 |
| `community_network_modularity.csv` | 05 | Walktrap modularity per network in the algorithm comparison |
| `badger_prescribed_groups.csv` | 05 | Badger field-defined groups and their predicted classes — **not** part of the algorithm comparison (Supp. Table S12) |
| `community_method_failures.csv` | 05 | (graph, method) pairs excluded because the algorithm did not converge |
| `table_s6_representative.csv` | 05 | Representative subset behind Supp. Table S6 |
| `threshold_sensitivity_{raw,by_graph,wide}.csv` | 06 | Class stability under edge pruning |
| `threshold_sensitivity_by_class.csv` | 06 | Headline invariance by class (Table 3; by taxon, Supp. Table S1) |
| `table_s14_ants.csv` | 08 | Ant networks and groups (Supp. Table S14) |
| `metadata_sensitivity_{overall,classwise,auc}.csv` | supporting | Metadata-model sensitivity checks |

### Where each table appears in the manuscript

| Manuscript table | Source |
|------------------|--------|
| Table 2 — first-stage classification performance | Adapted from Appaw et al. (2025) |
| Table 3 — threshold sensitivity by class | `threshold_sensitivity_by_class.csv` (06) |
| Table 4 — representative networks and their groups | `merged_community_predictions.csv` (05) |
| Tables 5 & 6 — metadata models vs majority baseline | `metrics_meta1.csv`, `metrics_meta2.csv` (03) |
| Table 7 — grouped cross-validation | `pipeline_overall_performance2.csv` (04) |
| Supp. S1 — threshold sensitivity by taxon | `threshold_sensitivity_by_graph.csv` (06) |
| Supp. S2–S5 — category sample sizes | `sample_sizes_by_*.csv` (02) |
| Supp. S6 — cross-algorithm group classification | `table_s6_representative.csv` (05) |
| Supp. S7 — metadata confusion matrices | `classwise_meta1.csv`, `classwise_meta2.csv` (03) |
| Supp. S8 — full grouped-CV breakdown | `pipeline_overall_performance2.csv` (04) |
| Supp. S12 — badger | `badger_prescribed_groups.csv` (05) |
| Supp. S13 — chimpanzee | `merged_community_predictions.csv` (05) |
| Supp. S14 — ant | `table_s14_ants.csv` (08) |

The single giraffe network is reported in a main-text table rather than in the
supplement; its rows are in `merged_community_predictions.csv` with the others.

**Two sources of ant rows.** `merged_community_predictions.csv` holds the ant
networks that enter the four-algorithm comparison (the six aggregate colonies
plus colony-1 day05 and day20). `table_s14_ants.csv` holds the full ant table
(the six aggregates plus days 05, 20 and 35 for every colony, walktrap only).
Supplementary Table S14 is built from the latter.

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

To re-export the portable boosters, run `recover_booster()` **once** in an
environment with the *older* XGBoost that wrote the original files (the calls are
at the top of `03_metadata_model_evaluation.R`, commented out).

Every script that needs the first-stage classifier — 02, 05, 06, 08 and
`sbm-sim-different-communities.R` — goes through `load_main_model()` /
`predict_main_model()`.

**Outputs generated before the recovery are not comparable with current ones.**
Any CSV in `outputs/` predating it was produced by a different model route and by
`01_functions.R` before the igraph 2.0 return-type fix in `classify_subgraphs()`.
Re-run the affected step rather than reusing an old file.

---

## Notes on three analyses

**Minimum group size (step 05).** Groups smaller than **8 nodes** are not
classified (`MIN_GROUP_SIZE` in `00_setup.R`). Field-defined badger social groups
span 2–9 individuals, so a threshold of 8 keeps the largest genuine groups and
excludes groups too small to yield stable structural features. The floor retains
62% of all detected groups. Networks with no qualifying group are reported as
such rather than being force-classified.

**Scope of the community-detection comparison (step 05) — badgers are not in it.**
The comparison across Walktrap, Louvain, leading-eigenvector and Infomap runs on
the **chimpanzee, giraffe and ant networks only**. Badger networks are excluded
by design: their groups are field-defined (Silk, Weber et al. 2018) rather than
detected, so running detection algorithms on them cannot test whether *detected*
groups are stable. Badgers are classified separately from their prescribed
groups, written to `badger_prescribed_groups.csv`, and their rows carry
`method = "prescribed"` in `merged_community_predictions.csv`. They contribute
nothing to `community_method_stability.csv`, which is the cross-algorithm result.

In the code this is the split between `DETECTION_NETWORKS` and
`BADGER_NETWORKS` at the top of `05_group_level_community_detection.R`. If you
add networks, put them in whichever list matches how their groups were obtained —
adding a prescribed-group network to `DETECTION_NETWORKS` would silently
contaminate the stability figures.

The ant set contains the six aggregate colony networks **and** the two daily
colony-1 networks with the highest modularity. The aggregate colony networks have
negligible modularity (< 0.01), so any partition of them is an arbitrary division
of an essentially unstructured graph; the daily networks carry real community
structure and are the meaningful test for ants. The two daily networks were
selected by **modularity**, not by predicted class — selecting on the outcome
would be circular. `community_network_modularity.csv` reports the modularity of
every network in the comparison, and `community_method_stability.csv` reports the
geometric percentages both for all networks and restricted to those with
modularity ≥ 0.01. Note that the modularity floor also removes nine chimpanzee
networks, so the restricted subset is "networks with detectable community
structure", not an ant-specific correction.

Infomap returns a single community spanning the whole network for every network
in the comparison, so its percentage is a network-level rather than a group-level
figure and is not directly comparable with the other three methods.

**Threshold sensitivity (step 06).** Most ASNR networks carry edge weights, but
the meaning of those weights varies among studies — interaction frequency,
duration or intensity — so they are not directly comparable across species and
study systems. Following Appaw et al. (2025) we therefore treat every network as
static and binarised, and restrict the analysis to the largest connected
component. **Binarisation is a modelling choice, not a property of the data**;
the first-stage classifier and its simulated training networks are defined on
unweighted graphs, so a weighted re-analysis would require the classifier to be
retrained rather than simply re-applied.

Because the pipeline runs on binarised networks, thresholding on observed
interaction weights is not available within the present model set. Instead we
remove edges by a structural connection-strength proxy (neighbourhood overlap /
Jaccard of shared neighbours) and re-classify at each step. This is a
**topological robustness check** — it is not a re-analysis of contact frequency
or duration, and should not be described as one.

---

## Environment

- Key packages: `igraph`, `tidymodels` (`parsnip`, `recipes`, `workflows`, `tune`,
  `finetune`, `rsample`, `yardstick`), `xgboost`, `shapviz`, `hstats`,
  `missForest`, `themis`, `janitor`, `here`, `tidyverse`

Exact versions are pinned in `renv.lock` (242 packages); `renv::restore()`
reproduces them.

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
