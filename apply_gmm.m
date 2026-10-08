%% EXAMPLE: APPLY THE TRAINED GMM TO NEW MEALS
% Applies the trained model (gmm_model.mat) to example.csv using detectSAmeals function.
% To use your own data, replace "example.csv" with the path of your CSV
% (columns RaAUC120min, TGR50, d, not normalized; optional ID).
% Requires: Statistics and Machine Learning Toolbox
clear all, close all, clc


results = detectSAmeals("example/example.csv", ...
    OutputFile = "example_assignments.csv", ...
    PlotPosterior = true, ...
    PlotMixtures = true);

results
