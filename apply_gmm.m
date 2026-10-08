%% APPLY THE TRAINED GMM TO NEW MEALS
% Assigns new meals to cluster FA or SA using the GMM trained by
% analysis_and_plot_generation.m (gmm_model.mat). New data are normalized
% with the mean and SD of the training set.
% Requires: Statistics and Machine Learning Toolbox
% Input: CSV with one row per meal and (at least) the columns used for the
% clustering: RaAUC120min, TGR50, d (not normalized). An "ID" column is optional.
% Output: table with posterior probabilities and labels, saved as CSV
clear all, close all, clc

addpath("functions\")

%% Settings
inputFile  = 'example.csv';
modelFile  = 'gmm_model.mat';
outputFile = 'assignment_results.csv';

%% Loading model and new data
load(['model\' modelFile],'gmm_model')
% gmm_model is a structure with fields:
% - gm_dist: the trained Gaussian Mixture Model (gmdistribution object)
% - featNames: cell array with feature names that gm_dist requires in input
% - mu: training-set mean to z-score the new observations 
% - sigma: training-set standard deviation to z-score the new observations
% - idxFA & idxSA: indeces of clusters FA and SA
% - threshold: uncertain posterior lower and upper bound

newTable = readtable(['example\' inputFile]);

%% Normalizing with the training mean and SD
X = table2array(newTable(:,gmm_model.featNames));
Xz = (X - gmm_model.mu)./gmm_model.sigma;

%% Posterior probability
post = posterior(gmm_model.gm_dist,Xz);

% defining cluster labels
idxFA = gmm_model.idxFA; idxSA = gmm_model.idxSA;
postFA = post(:,idxFA);
postSA = post(:,idxSA);

%% Labels (maximum posterior) and uncertain meals
label = strings(size(Xz,1),1);
isSA = false(size(Xz,1),1);
[~,idxk_new] = max(post,[],2);
isSA = idxk_new == idxSA;
label(1:size(X,1)) = "FA";
label(isSA) = "SA";

uncertain = postSA >= gmm_model.threshold(1) & postSA <= gmm_model.threshold(2);

results = table(postFA,postSA,label,uncertain, ...
    'VariableNames',{'Posterior_FA','Posterior_SA','Label','Uncertain'});
if ismember('ID',newTable.Properties.VariableNames)
    results = [newTable(:,'ID') results];
end

results
writetable(results,['example\' outputFile])

disp([num2str(sum(label=="SA")) ' meals assigned to SA, ' num2str(sum(label=="FA")) ...
    ' to FA'])

%% Plot: posterior probability in the (normalized) feature space
plotPosterior(gmm_model,Xz,label)
%% Plot: observation assignment overlayed on the Gaussian Mixtures projected onto each 2D-feature pair
plotMixtures(gmm_model,Xz)
