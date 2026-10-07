%% HYPERPARAMETER TUNING: FEATURES, CLUSTERING ALGORITHM AND NUMBER OF CLUSTERS
% Beam-search feature selection on the training set. For every candidate
% feature set, k-means, GMM and H-clust are evaluated (number of clusters
% chosen by the silhouette) and the best algorithm is retained. The selected
% set is then clustered with each algorithm and the agreement between
% the algorithms is computed (NMI).
% Requires: Statistics and Machine Learning Toolbox
% Input (working directory): PARAMETERS_notNormalized.csv, functions/ folder
clear all, close all, clc
warning('off')
addpath("functions")

paramsTable = readtable("PARAMETERS_notNormalized.csv");

% Selecting DDF
% paramsTable = paramsTable(:,[19:33]); % DDF only

% Or Selecting MBF+DDF (excluding demographics and correlated features)
paramsTable = paramsTable(:,[6 8 10:33]); % MBF + DDF

fullTABLE = readtable("PARAMETERS_notNormalized.csv");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       CREATE TEST SET                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1)
ID = table2array(fullTABLE(:,"ID"));
uniqueID = unique(ID);
cv = cvpartition(length(uniqueID),'HoldOut',0.3);
trainID = uniqueID(training(cv));
testID = uniqueID(test(cv));

idxTraining = ismember(ID, trainID);
idxTest  = ismember(ID, testID);

paramsTable = paramsTable(idxTraining,:);
fullTABLE = fullTABLE(idxTraining,:);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_notNormalized = table2array(paramsTable);
params = zscore(params_notNormalized);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        BEAM SEARCH FEATURE SELECTION (GMM / k-means / H-clust)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

Nfeatures = size(paramsTable,2);
parNames = paramsTable.Properties.VariableNames;

%% Search settings
beamWidth         = 5;      % number of feature sets kept at every step (1 -> plain forward selection)
maxFeatures       = 7;      % limit on maximum number of features 
minRelImprovement = 0.005;  % significant improvement only if current best score > best*(1+0.005)
patience          = 1;      % non-improving steps tolerated before stopping
scoreCeiling      = 0.9999;  % stop if the score is almost perfect

%% Initialization
beam = zeros(1,0);          % each row is one feature set (indices)
best_score = 0;
best_features_by_score = []; % indeces of best feature set
best_alg = ''; best_bestk = NaN; 
best_beam_sets = []; best_beam_scores = [];   % book-keeping the beam at the best score 
stall = 0;                  % consecutive steps without improvement

% auxiliary variables for beam search hyperparameters decision process
path_size = []; path_features = {}; path_score = []; path_alg = {}; path_k = [];

for step = 1:maxFeatures

    % --- expand every set in the beam with every missing feature ---------
    C = zeros(0,step);
    for b = 1:size(beam,1)
        base = beam(b,:);
        for f = setdiff(1:Nfeatures, base)
            C(end+1,:) = sort([base f]); %#ok<SAGROW>
        end
    end
    C = unique(C,'rows');   % delete possible duplicate combinations

    % --- evaluate all unique  feature set candidates ---------------------
    nCand = size(C,1);
    cand_score = zeros(nCand,1);
    cand_alg   = cell(nCand,1);
    cand_k     = zeros(nCand,1);

    for c = 1:nCand  
        [cand_score(c), cand_alg{c}, cand_k(c)] = ...
            evaluateSet(C(c,:), params, params_notNormalized, fullTABLE);
    end

    % --- keep the best sets ----------------------------------------------
    [cand_score, order] = sort(cand_score,'descend');
    C = C(order,:); cand_alg = cand_alg(order); cand_k = cand_k(order);
    nKeep = min(beamWidth, nCand);
    beam = C(1:nKeep,:);

    % storing best score and best set at current step
    step_score = cand_score(1);
    step_set   = beam(1,:);

    % book-keeping decision process
    path_size(end+1)     = step;                              %#ok<SAGROW>
    path_features{end+1} = strjoin(parNames(step_set), ', '); %#ok<SAGROW>
    path_score(end+1)    = step_score;                        %#ok<SAGROW>
    path_alg{end+1}      = cand_alg{1};                       %#ok<SAGROW>
    path_k(end+1)        = cand_k(1);                         %#ok<SAGROW>

    % --- improvement check against the current best set -------------
    if step_score > best_score * (1 + minRelImprovement)
        best_score = step_score;
        best_features_by_score = step_set;
        best_alg = cand_alg{1};
        best_bestk = cand_k(1);
        best_beam_sets = beam;
        best_beam_scores = cand_score(1:nKeep);
        stall = 0;
        status = 'improved';
    else
        stall = stall + 1;
        status = sprintf('no improvement (%d/%d)', stall, patience);
    end

    fprintf('\nStep %d | top score %.4f | %s\n', step, step_score, status);
    for b = 1:nKeep
        fprintf('   #%d  %.4f  [%s, k=%d]  %s\n', b, cand_score(b), cand_alg{b}, cand_k(b), ...
            strjoin(parNames(beam(b,:)), ', '));
    end

    % --- stopping rules --------------------------------------------------
    if best_score + best_score*minRelImprovement > 1
        disp('Next best score required to be > 1. Stopping search...')
        break
    end

    if best_score >= scoreCeiling
        disp('Stop: score ceiling reached')
        break
    end
    if stall > patience
        disp('Stop: patience exhausted (no relative improvement > threshold)')
        break
    end
end

% --- best set  -----------------------------------------------------------
disp('=========== SELECTED FEATURES ===========')
disp(parNames(best_features_by_score))
disp(['Score: ' num2str(best_score)])
disp(['#' num2str(best_bestk) ' cluster solution - with ' best_alg])
pathTable = table(path_size', path_features', path_score', path_alg', path_k', ...
    'VariableNames',{'Size','TopSet','Score','BestAlgorithm','K'})

disp('Ranking of best sets found:')
for b = 1:size(best_beam_sets,1)
    fprintf('   #%d  %.4f  %s\n', b, best_beam_scores(b), strjoin(parNames(best_beam_sets(b,:)), ', '));
end


%%
parset_notNormalized = params_notNormalized(:,best_features_by_score);
parset = params(:,best_features_by_score);

disp('Press any button to analyze CV with the best feature set')
pause()
clc

%% ====================== kmeans ======================= %%
% k-means
disp('===================== k-MEANS =======================')
rng(1)
[~, bestk,meanS,sdS] = kmeans_bestk(parset);
rng(1)
idxk_kmeans = kmeans(parset,bestk,"Replicates",5,"Distance","sqeuclidean",'Display','off',MaxIter=1e4);
parset_prediction = cv_kmeans_with_SD(fullTABLE,{'Fat_amt';'Prot_amt'},bestk,parset_notNormalized);
counts = histcounts(idxk_kmeans,1:bestk+1);
p_kmeans = counts/sum(counts);
entropy = -sum(p_kmeans.*log2(p_kmeans))/log2(bestk);
k_means_score = mean(table2array(parset_prediction([1 2],2)))*entropy;
mc_agreement = mean(table2array(parset_prediction([1 2], 2)));

fprintf('\n=== Clustering Summary ===\n');
fprintf('  %-18s %8.4f\n', 'MC-rank agreement:', mc_agreement);
fprintf('  %-18s %8.4f\n', 'Entropy:',           entropy);
fprintf('  %-18s %8.4f\n', 'Score:',         k_means_score);
fprintf('==========================\n');
fprintf('Moving to GMM...\n\n');
pause()
clc


%% ======================= GMM ======================== %%
disp('===================== GMM =======================')
[bestk,meanS,sdS] = gmm_bestk(parset,6);
parset_prediction = cv_gm_with_SD(fullTABLE,{'Fat_amt';'Prot_amt'},bestk,parset_notNormalized);
curr_gm_dist = fitgmdist(parset,bestk, ...
    'RegularizationValue',0.01,'CovarianceType','diagonal','SharedCovariance',true);
% hard-assignment based on the gaussian component with the maximum posterior
post = posterior(curr_gm_dist,parset);
[~,idxk_gmm] = max(post,[],2);
p_gmm = curr_gm_dist.ComponentProportion;
entropy = -sum(p_gmm.*log2(p_gmm))/log2(bestk);
gmm_score = mean(table2array(parset_prediction([1 2],2)))*entropy;
mc_agreement = mean(table2array(parset_prediction([1 2], 2)));

fprintf('\n=== Clustering Summary ===\n');
fprintf('  %-18s %8.4f\n', 'MC-rank agreement:', mc_agreement);
fprintf('  %-18s %8.4f\n', 'Entropy:',           entropy);
fprintf('  %-18s %8.4f\n', 'Score:',         gmm_score);
fprintf('==========================\n');
fprintf('Moving to H-clust...\n\n');
pause()
clc

%% ======================= H-clust ======================== %%
disp('===================== H-clust =======================')
[~, bestk, meanS, sdS] = hclust_bestk(parset);
h = linkage(parset,"ward");
idxk_h = cluster(h,"maxclust",bestk);
parset_prediction = cv_hclust_with_SD(fullTABLE,{'Fat_amt';'Prot_amt'},bestk,parset_notNormalized);
counts = histcounts(idxk_h,1:bestk+1);
p_h = counts/sum(counts);
entropy = -sum(p_h.*log2(p_h))/log2(bestk);
h_score = mean(table2array(parset_prediction([1 2],2)))*entropy;
mc_agreement = mean(table2array(parset_prediction([1 2], 2)));

fprintf('\n=== Clustering Summary ===\n');
fprintf('  %-18s %8.4f\n', 'MC-rank agreement:', mc_agreement);
fprintf('  %-18s %8.4f\n', 'Entropy:',           entropy);
fprintf('  %-18s %8.4f\n', 'Score:',         h_score);
fprintf('==========================\n');
fprintf('Moving to H-clust...\n\n');
pause()
clc

disp('=== Normalized Mutual Information ===')
disp(['NMI between kmeans and h-clust: ' num2str(computeNMI(idxk_kmeans,p_kmeans,idxk_h,p_h))])
disp(['NMI between kmeans and gmm: ' num2str(computeNMI(idxk_kmeans,p_kmeans,idxk_gmm,p_gmm))])
disp(['NMI between h-clust and gmm: ' num2str(computeNMI(idxk_h,p_h,idxk_gmm,p_gmm))])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          LOCAL FUNCTIONS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [score, alg, bestk_out] = evaluateSet(idx, params, params_notNormalized, fullTABLE)
    % Score of a feature set = best (over k-means, GMM, H-clust) of
    % CV accuracy x normalized cluster-size entropy.
    str = {'Fat_amt','Prot_amt'};
    parset = params(:,idx);
    parset_notNormalized = params_notNormalized(:,idx);

    % ---- k-means ----
    rng(1)
    [~, bestk_km] = kmeans_bestk(parset);
    rng(1)
    idxk = kmeans(parset,bestk_km,"Replicates",5,"Distance","sqeuclidean",'Display','off',MaxIter=1e4);
    pred = cv_kmeans(fullTABLE,str,bestk_km,parset_notNormalized);
    counts = histcounts(idxk,1:bestk_km+1);
    p = counts/sum(counts);
    ent = -sum(p.*log2(p))/log2(bestk_km);
    s_km = mean(pred(1:2))*ent;

    % ---- GMM ----
    bestk_gm = gmm_bestk(parset,6);
    pred = cv_gm(fullTABLE,str,bestk_gm,parset_notNormalized);
    gm = fitgmdist(parset,bestk_gm, ...
        'RegularizationValue',0.01,'CovarianceType','diagonal','SharedCovariance',true);
    p = gm.ComponentProportion;
    ent = -sum(p.*log2(p))/log2(bestk_gm);
    s_gm = mean(pred(1:2))*ent;

    % ---- H-clust ----
    [~, bestk_h] = hclust_bestk(parset);
    h = linkage(parset,"ward");
    idxk = cluster(h,"maxclust",bestk_h);
    pred = cv_hclust(fullTABLE,str,bestk_h,parset_notNormalized);
    counts = histcounts(idxk,1:bestk_h+1);
    p = counts/sum(counts);
    ent = -sum(p.*log2(p))/log2(bestk_h);
    s_h = mean(pred(1:2))*ent;

    scores = [s_km, s_gm, s_h];
    [score, a] = max(scores);
    algs = {'kmeans','gmm','hclust'};
    ks   = [bestk_km, bestk_gm, bestk_h];
    alg = algs{a};
    bestk_out = ks(a);
end


function NMI = computeNMI(assignments_C1, proportions_C1, assignments_C2, proportions_C2)
    % Normalized mutual information (geometric-mean normalization) between
    % two clusterings, given their assignments and cluster proportions.
    I_C1_C2 = 0;
    % entropy of C1 and C2
    entropy_C1 = -sum(proportions_C1.*log2(proportions_C1));
    entropy_C2 = -sum(proportions_C2.*log2(proportions_C2));

    n_groups_C1 = length(unique(assignments_C1));
    n_groups_C2 = length(unique(assignments_C2));

    tot_observations = length(assignments_C1);

    % joint proportions (pxy)
    joint_counts = accumarray([assignments_C1 assignments_C2], 1, [n_groups_C1 n_groups_C2]);
    joint_proportions = joint_counts/tot_observations;

    % outer proportions (px_py)
    outer_proportions = proportions_C1'*proportions_C2;
    
    for ii = 1:size(joint_proportions,1)
        for jj = 1:size(joint_proportions,2)
            if joint_proportions(ii,jj) > 0
                % mutual information
                 I_C1_C2 = I_C1_C2 + joint_proportions(ii,jj)*log2(joint_proportions(ii,jj) / outer_proportions(ii,jj));
            end
        end
    end
    % normalize
    NMI = I_C1_C2 / sqrt(entropy_C2*entropy_C1);
end