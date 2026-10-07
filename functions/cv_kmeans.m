function [prediction_performance] = cv_kmeans(fullTABLE,str,bestk,params)
% CV_KMEANS performs 5-fold subject-level cross-validation of the k-means
% clustering using k = bestk and the parameters in params.
% A reference clustering is computed on all the data; in each fold k-means
% is fitted on the training subjects only and the test observations are
% assigned to the closest training centroid.
% INPUT:
% fullTABLE: table with the (not normalized) data, including the "ID" column
% str: cell array with the names of the macronutrient variables of interest
% bestk: number of clusters
% params: matrix of the (not normalized) parameters used for clustering
% OUTPUT:
% prediction_performance: #macronutrients x 1 fraction of test observations
% assigned to the same (ordered) cluster as in the reference clustering

%% creating reference clustering
rng(1)
params_normalized = zscore(params);
idxRef = kmeans(params_normalized,bestk,"Replicates",5,"Distance","sqeuclidean",MaxIter=1e4);

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

    %% K-means on training data
    [idxCV,trainCentroids] = kmeans(currTrainSet,bestk,"Replicates",5,"Distance","sqeuclidean",MaxIter=1e4);

    %% median of the current clustering (using ONLY training data)
    currTABLE = normalize(fullTABLE(currTrainIdx,:),"zscore");
    currMedian = clustersMedian(currTABLE,str,bestk,idxCV);
    [~,orderCV] = sort(currMedian,2);
    
    %% assign test data to the closest centroid
    currCentroidDist = pdist2(currTestSet, trainCentroids, 'squaredeuclidean');
    [~, assignedIdx] = min(currCentroidDist, [], 2);

    newAssignedIdx = zeros(length(assignedIdx),1);
    newIdxRef = zeros(length(idxRef),1);

    %% ordering cluster indices with the sorted medians
    currSum = zeros(length(str),1);
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
