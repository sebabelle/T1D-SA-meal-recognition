function prediction_performance = cv_gm(fullTABLE,str,bestk,params)
% CV_GM performs 5-fold subject-level cross-validation of the GMM clustering
% (diagonal, shared covariance) using bestk components and the parameters in
% params.
% A reference clustering is computed on all the data; in each fold the GMM
% is fitted on the training subjects only and the test observations are
% assigned to the component with the maximum posterior.
% INPUT:
% fullTABLE: table with the (not normalized) data, including the "ID" column
% str: cell array with the names of the macronutrient variables of interest
% bestk: number of components/clusters
% params: matrix of the (not normalized) parameters used for clustering
% OUTPUT:
% prediction_performance: #macronutrients x 1 fraction of test observations
% assigned to the same (ordered) cluster as in the reference clustering

%% creating reference clustering
rng(1)
params_normalized = zscore(params);
gm_dist_ref = fitgmdist(params_normalized,bestk, ...
    'RegularizationValue',0.01,'CovarianceType','diagonal','SharedCovariance',true);
idxRef = cluster(gm_dist_ref,params_normalized);

%% matching reference and CV clusters
% Cluster labels are arbitrary and can differ between the reference and the
% CV clusterings, so clusters are matched by ordering them according to the
% median of each macronutrient (lowest median -> 1, ..., highest -> bestk)
refMedian = clustersMedian(normalize(fullTABLE,"zscore"),str,bestk,idxRef);
[~,orderRef] = sort(refMedian,2);

%% creating 5-fold subject-level CV
ID = table2array(fullTABLE(:,"ID"));
uniqueID = unique(ID);
cv = cvpartition(length(uniqueID),'Kfold',5);

% sum of the correct assignments over all folds
tmp_sum = zeros(length(str),1);

for fold = 1:5 % for each fold
    % training and test indices of the current fold
    currTrainID = uniqueID(training(cv,fold));
    currTestID = uniqueID(test(cv,fold));

    currTrainIdx = ismember(ID,currTrainID);
    currTestIdx  = ismember(ID,currTestID);

    % normalizing the test set with the training mean and SD
    [currTrainSet, curr_mu_train, curr_sigma_train] = zscore(params(currTrainIdx,:));
    currTestSet = (params(currTestIdx,:) - curr_mu_train)./curr_sigma_train;

    %% GMM clustering on training set
    gm_dist_train = fitgmdist(currTrainSet,bestk, ...
        'RegularizationValue',0.01,'CovarianceType','diagonal','SharedCovariance',true);
    idxCV = cluster(gm_dist_train,currTrainSet);

    %% median of the current clustering (using ONLY training data)
    currTABLE = normalize(fullTABLE(currTrainIdx,:),"zscore");
    currMedian = clustersMedian(currTABLE,str,bestk,idxCV);
    [~,orderCV] = sort(currMedian,2);
    
    %% assign test data to the maximum posterior
    currPosterior = posterior(gm_dist_train,currTestSet);
    [~,assignedIdx] = max(currPosterior,[],2);

    %% ordering cluster indices with the sorted medians
    currSum = zeros(length(str),1);
    newIdxRef = zeros(length(idxRef),1);
    newAssignedIdx = zeros(length(assignedIdx),1);
    for macronutrient = 1:length(str) % for each macronutrient
        for kk = 1:bestk % for each cluster
            newIdxRef(idxRef == orderRef(macronutrient,kk)) = kk;
            newAssignedIdx(assignedIdx == orderCV(macronutrient,kk)) = kk;
        end
        % correct predictions of the current fold
        currSum(macronutrient,:) = sum(newIdxRef(currTestIdx) == newAssignedIdx);
        tmp_sum(macronutrient,:) = tmp_sum(macronutrient,:) + currSum(macronutrient,:);
    end
end

prediction_performance = tmp_sum./length(idxRef);

end
