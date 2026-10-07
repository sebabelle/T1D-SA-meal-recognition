# T1D-SA-meal-recognition
Code supporting the manuscript:

> **A Physiology-Informed Approach for Automated Recognition of Slow-Absorption Meals Under Free-Living Conditions in Type 1 Diabetes**, Bellese et al.

## Requirements
- MATLAB R2025b or newer
- Statistics and Machine Learning Toolbox

## Ready-to-use model and example of usage

This repository provides the **trained Gaussian mixture model** (`gmm_model.mat`) and an **example input file**
(`example.csv`).

**Quick start**

1. Open MATLAB in the repository folder.
2. Run `apply_gmm.m`. By default it applies `gmm_model.mat` to `example.csv`.
3. The script displays the results, saves them to `example_assignments.csv` and plots the posterior probabilities
   and the fitted mixture components together with the meals.

**Using your own data**: prepare a CSV like `example.csv`, with one row per meal and the columns

| Column | Description |
|---|---|
| `RaAUC120min`, `TGR50`, `d` | Clustering features, **not normalized**, computed as in the manuscript |
| `ID` | Subject identifier |

then set `inputFile` and `outputFile` at the top of `apply_gmm.m`.

**Output**: one row per meal with `Posterior_FA`, `Posterior_SA`, `Label` (`FA` or `SA`, maximum posterior) and
`Uncertain` (posterior of SA between 0.3 and 0.7).
Note: data are normalized internally with the mean and SD of the training set used to developed the model.

**About the model**: The model is a GMM with two components with diagonal shared covariance, and regularization 0.01. `gmm_model.mat` contains the fitted `gmdistribution` MATLAB object, the feature names, the training
mean/SD used for normalization, the cluster indices and the uncertainty thresholds. It contains no participant data.
It may be regenerated with `analysis_and_plot_generation.m`.



## Repository structure

```
hyperparameter_tuning.m             1. selects features, clustering algorithm and number of clusters
analysis_and_plot_generation.m      2. trains and tests the Gaussian mixture model, performs analyses and generates the plots presented in the manuscript
apply_gmm.m                         3. example of use of the trained Gaussian mixture model to detect SA meals and visualization examples
functions/                          helper functions
model/                              contains the trained and ready-to-use Gaussian mixture model
example/                            contains a small synthetic dataset ready to be tested on apply_gmm.m
 
README.md
```

1. **`hyperparameter_tuning.m`**: Performs subject-level hold-out split (70% training / 30% test), then a
   beam-search forward feature selection on the training set based on a unsupervised subject-wise 5-fold CV procedure. Each candidate feature set is scored as the best, over
   k-means, GMM and H-clust, of `mean(CV accuracy for Fat and Protein) x normalized entropy of the cluster proportions`;
   the number of clusters is chosen by maximization of the mean silhouette (k = 2-6). The beam search is regulated by `beamWidth`, `maxFeatures`,`minRelImprovement`,`patience`,`scoreCeiling`.      
   To reproduce the result reported in the manuscript, use `beamWidth=5`,
   `maxFeatures=7`, `minRelImprovement=0.005`, `patience=1`, `scoreCeiling=0.9999`.
   The agreement between algorithms is computed (NMI). The test set is not used.
3. **`analysis_and_plot_generation.m`**: Using the selected features, fits the GMM  on the training subjects, computes the CV accuracy,
   entropy, Calinski-Harabasz and Davies-Bouldin indices, Wilcoxon rank sum tests and the soft-clustering
   analyses, assigns the test subjects through the posterior probability and generates the figures reported in the manuscript
   (3D scatter of the posterior, boxcharts, median CGM profiles). Saves `gmm_model.mat` (the Gaussian Mixture model), `assignments.mat` (meal labels assigned by the GMM),
   `posterior.mat` (posterior probability of cluster membership).
4. **`apply_gmm.m`**: Loads `gmm_model.mat` and a CSV of new meals to classify, normalizes them with the training mean and SD,
   and returns the posterior probabilities, the label (`FA` / `SA`) and an `Uncertain` flag (posterior of `SA` / `FA` within
   0.3-0.7). Results are saved as CSV.


 ### Functions

| Function | Description |
|---|---|
| `kmeans_bestk`, `hclust_bestk`, `gmm_bestk` | Select the number of clusters (k = 2-6) maximizing the mean silhouette |
| `cv_kmeans`, `cv_hclust`, `cv_gm` | Perform the 5-fold subject-level cross-validation|
| `cv_kmeans_with_SD`, `cv_hclust_with_SD`, `cv_gm_with_SD` | Perform the 5-fold subject-level cross-validation, with per-fold output and SD across folds |
| `clustersMedian` | Computes the median of each macronutrient in each cluster (used to match cluster labels across clusterings) |
| `CH_index`, `DB_index` | Calinski-Harabasz and Davies-Bouldin indices (expect z-scored data) |
| `mahalanobisDistance` | Per-feature Mahalanobis distance between the means of two GMM components |
| `MacronutrientsSumRankTest` | Wilcoxon rank sum test between two clusters for Fat and Protein|
| `extractingCGMtracks`, `prctiles_curve` | Extract and plots (spaghetti plot) CGM profiles of a group of meals and returns median and 25th-75th percentile band |
| `plotPosterior`, `plotMixtures` | Plot the estimated posterior and the observation overlayed on the Gaussian Mixture components |

## Input data

For `hyperparameter_tuning.m` and `analysis_and_plot_generation.m` the following is needed:

- `PARAMETERS_notNormalized.csv`: one row per meal, with an `ID` column (subject identifier), the clustering
  parameters (including, at least, `RaAUC120min`, `TGR50`, `d`) and the macronutrient variables `Fat_amt` and `Prot_amt`,
   expressed in grams.
- `Dataset/`: one `.mat` file per meal, each with a struct `data` containing `t_CGM` (time grid) and `CGM`. The files must be in the same order as the rows of the CSV.

## Reproducibility

Random seeds are fixed with `rng(1)` for reproducibility. Results may differ slightly across MATLAB releases or toolbox versions.
