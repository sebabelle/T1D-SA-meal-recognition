# T1D-SA-meal-recognition
Code supporting the manuscript:

> **A Physiology-Informed Approach for Automated Recognition of Slow-Absorption Meals Under Free-Living Conditions in Type 1 Diabetes**, Bellese et al.

## Requirements
- MATLAB R2025b or newer
- Statistics and Machine Learning Toolbox

For more details see [Reproducibility Notes](https://github.com/sebabelle/T1D-SA-meal-recognition/blob/main/README.md#reproducibility).

## Introduction 

This repository provides the **trained Gaussian mixture model** (`gmm_model.mat`) inside the folder `model/` and the **MATLAB function** `detectSAmeals()` that allows to apply the GMM model on new meals. The functions assigns them to a class: slow-absorption (`SA`) or fast-absorption (`FA`) using three **model-based** features related to gastro-intestinal absorption.

**About the model**: The model is a GMM with two components with diagonal shared covariance, and regularization 0.01. `gmm_model.mat` contains the fitted `gmdistribution` MATLAB object, the feature names, the training
mean/SD used for normalization, the cluster indices and the uncertainty thresholds. It contains no training data.
It may be regenerated with `analysis_and_plot_generation.m`.

Part of this repository is dedicated to report the code used for **training and testing** of the **Gaussian mixture model** and **analysis** reported in the manuscript.

For a full description of the repository and its files see [Repository Structure](https://github.com/sebabelle/T1D-SA-meal-recognition/blob/main/README.md#repository-structure).

## How to use the model and example of usage


An example of **input file** (`example.csv`) and  **output file** (`assignment_results.csv`) is inside the folder `example/`.

### Quick start with the example:

1. Open MATLAB in the repository folder.
2. Run `apply_gmm.m`. By default it applies `gmm_model.mat` to `example.csv` through the function `detectSAmeals`.
3. The script displays the results, saves them to `example_assignments.csv` inside the folder `example\` and plots the posterior probabilities
   and the fitted mixture components together with the meals.

### Using your own data: 
1. Prepare a CSV file like `example.csv` (e.g.,`my_meals.csv`) , with one row per meal and the columns `RaAUC120min`, `TGR50`, `d`, and, optionally `ID`, as reported below. Alternatively a MATLAB table with the same structure can be used.

| Column | Description |
|---|---|
| `RaAUC120min`, `TGR50`, `d` | Clustering features, **not normalized**, computed as in the manuscript |
| `ID` | Subject identifier (Optional) |

| ID | RaAUC120min | TGR50 | d |
| ---: | ---: | ---: | ---: |
| 1 | 0.31 | 130.1 | 0.24 |
| 1 | 0.48 | 95.2 | 0.26 |
| 1 | 0.46 | 112.1 | 0.23 |
| 3 | 0.62 | 70.5 | 0.05 |
| 3 | 0.74 | 55.9 | 0.35 |
| 4 | 0.54 | 70.1 | 0.22 |
| 5 | 0.66 | 48.0 | 0.001 |

2. Call the function `detectSAmeals()`, passing the CSV (or MATLAB table) containing 1 or more meals (e.g., `my_meals.csv`). Here are reported some examples:

```matlab
results = detectSAmeals("my_meals.csv"); % applies to "my_meals.csv" the GMM model, returns a table with assignment results
results = detectSAmeals("my_meals.csv",OutputFile="my_assignments.csv"); % additionally saves the results in "my_assignments.csv"
results = detectSAmeals("my_meals.csv",PlotPosterior=true,PlotMixtures=true); % additionally plot posterior and 2D feature space overlayed to trained mixtures
```
Note: Optional fields `PlotPosterior` and `PlotMixtures` require `PlotPosterior.m` and `PlotMixtures.m` to be inside `functions` folder.

Note: `detectSAmeals()` assumes by default that `gmm_model.mat` is under the folder `model/`. If the GMM model is in another location, you can specify it using `ModelFile="path/to/model/gmm_model.mat"`.

```matlab
results = detectSAmeals("my_meals.csv",ModelFile="path/to/model/gmm_model.mat"); % applies to "my_meals.csv" the GMM, specifying location of the trained GMM. Returns a table with assignment results
```

**Output**: `results` contains one row per meal with `Posterior_FA`, `Posterior_SA`, `Label` (`FA` or `SA`, maximum posterior) and
`Uncertain` (posterior of SA between 0.3 and 0.7). An example of output is reported below.
Note: data are normalized internally with the mean and SD of the training set used to developed the model.

| Column | Description |
|---|---|
| `Posterior_FA`, `Posterior_SA` | Estimated posterior of cluster assignment (`FA` and `SA`) |
| `Label` | Meal assigned label (`FA` or `SA`) |
| `Uncertaine` | Flag that return wheter posterior is uncertain (>0.3 & <0.7) |
| `ID` | Subject identifier (Optional) |

| ID | Posterior_FA | Posterior_SA | Label | Uncertain |
| ---: | ---: | ---: | :---: | ---: |
| 1 | 0.000333 | 0.999667 | SA | 0 |
| 1 | 0.314160 | 0.685840 | SA | 1 |
| 1 | 0.068700 | 0.931300 | SA | 0 |
| 3 | 0.998767 | 0.001233 | FA | 0 |
| 3 | 0.999784 | 0.000216 | FA | 0 |
| 4 | 0.966883 | 0.033117 | FA | 0 |
| 5 | 0.999963 | 0.000037 | FA | 0 |


## Repository structure

```
detectSAmeals.m                     1. function: applies the trained GMM to new meals (CSV or MATLAB table)
apply_gmm.m                         2. example of function call
hyperparameter_tuning.m             3. selects features, clustering algorithm and number of clusters
analysis_and_plot_generation.m      4. trains and tests the Gaussian mixture model, performs analyses and generates the plots presented in the manuscript

functions/                          helper functions
model/                              contains the trained and ready-to-use Gaussian mixture model
example/                            contains a small synthetic dataset ready to be tested on apply_gmm.m
 
README.md
```
1. **`detectSAmeals.m`**: Function that allows to apply the trained GMM to a new dataset. Loads `gmm_model.mat` and a CSV (or MATLAB table) of new meals, normalizes them with the training mean and SD,
   and returns the posterior probabilities, the label (`FA` / `SA`) and an `Uncertain` flag (posterior of SA within
   0.3-0.7). Optional fields: `PlotPosterior` and `PlotMixtures`plot the estimated posterior onto the 3D feature space and plots the fitted mixture components (1-2 SD contours) with the new meals. Results can be saved as CSV if specified.
2. **`apply_gmm.m`**: Reports an example of usage of the function `detectSAmeals.m`
3. **`hyperparameter_tuning.m`**: Performs subject-level hold-out split (70% training / 30% test), then a
   beam-search forward feature selection on the training set based on a unsupervised subject-wise 5-fold CV procedure. Each candidate feature set is scored as the best, over
   k-means, GMM and H-clust, of `mean(CV accuracy for Fat and Protein) x normalized entropy of the cluster proportions`;
   the number of clusters is chosen by maximization of the mean silhouette (k = 2-6). The beam search is regulated by `beamWidth`, `maxFeatures`,`minRelImprovement`,`patience`,`scoreCeiling`.      
   To reproduce the result reported in the manuscript, use `beamWidth=5`,
   `maxFeatures=7`, `minRelImprovement=0.005`, `patience=1`, `scoreCeiling=0.9999`.
   The agreement between algorithms is computed (NMI). The test set is not used.
4. **`analysis_and_plot_generation.m`**: Using the selected features, fits the GMM  on the training subjects, computes the CV accuracy,
   entropy, Calinski-Harabasz and Davies-Bouldin indices, Wilcoxon rank sum tests and the soft-clustering
   analyses, assigns the test subjects through the posterior probability and generates the figures reported in the manuscript
   (3D scatter of the posterior, boxcharts, median CGM profiles). Saves `gmm_model.mat` (the Gaussian Mixture model), `assignments.mat` (meal labels assigned by the GMM),
   `posterior.mat` (posterior probability of cluster membership).

To reproduce the results reported in the manuscript: run 3->4. To use the trained GMM model to classify new meals: run 1 only.


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
| `plotPosterior`, `plotMixtures` | Plot the estimated posterior and each observation overlayed on the Gaussian Mixture components |

## Input data required for model training and testing

For `hyperparameter_tuning.m` and `analysis_and_plot_generation.m` the following is needed:

- `PARAMETERS_notNormalized.csv`: one row per meal, with an `ID` column (subject identifier), the clustering
  parameters (including, at least, `RaAUC120min`, `TGR50`, `d`) and the macronutrient variables `Fat_amt` (Fat load) and `Prot_amt` (Protein load),
   expressed in grams.
- `Dataset/`: one `.mat` file per meal, each with a struct `data` containing `t_CGM` (time grid, [min]) and `CGM` (CGM samples, [mg/dl]). Must contain data within 1 hour before the mealtime up to 4 hours after. `t_CGM` is centered around mealtime (i.e., `t_CGM(mealtime)=0`) The files must be in the same order as the rows of the CSV.

## About Reproducibility

Random seeds are fixed with `rng(1)` for reproducibility across runs. 

Results may differ across MATLAB releases or toolbox versions.

Note: The results reported in the manuscript were obtained with MATLAB Version 25.2 R2025b Update 4. While the code may work with previous MATLAB releases, older versions may behave differently, which can arise errors or unexpected results.

