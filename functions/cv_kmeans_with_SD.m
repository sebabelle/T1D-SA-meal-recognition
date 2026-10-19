function [agreement] = cv_kmeans_with_SD(fullTABLE,str,bestk,params)
% CV_KMEANS_WITH_SD same as CV_KMEANS, but it also displays the per-fold
% results and returns the SD of the fold agreement.
% INPUT:
% fullTABLE: table with the (not normalized) data, including the "ID" column
% str: cell array with the names of the macronutrient variables of interest
% bestk: number of clusters
% params: matrix of the (not normalized) parameters used for clustering
% OUTPUT:
% agreement: table with, for each macronutrient, the mean
% agreement (fraction) and the SD of the fold agreement (%)

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

disp(' --- Performing Kfold Cross Validation --- ')

%% creating 5-fold subject-level CV
ID = table2array(fullTABLE(:,"ID"));
uniqueID = unique(ID);
cv = cvpartition(length(uniqueID),'Kfold',5);

% sum of the assignments in agreement with the reference clustering over all folds
tmp_sum = zeros(length(str),1);

% agreement of each fold
fold_agreement = zeros(5,length(str));

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
        % assignments in agreement with the reference clustering in the current fold
        currSum(macronutrient,:) = sum(newIdxRef(currTestIdx) == newAssignedIdx);
        tmp_sum(macronutrient,:) = tmp_sum(macronutrient,:) + currSum(macronutrient,:);
        fold_agreement(fold,macronutrient) = (currSum(macronutrient,:)/size(currTestSet,1))*100; 
    end

    %% displaying current fold agreement
    disp(['Analysing Fold n.' num2str(fold) '| Observations in agreement with the reference clustering: '])
    for mm = 1:length(str)
        disp([ str{mm} ': ' num2str(currSum(mm))  ' out of ' num2str(length(currTestSet)) ' | Fold Agreement ' num2str(fold_agreement(fold,mm)) '%'])
    end
    disp(' ------ ')
end

%% Creating agreement table
agreement = tmp_sum./length(idxRef);
std_dev_performance = std(fold_agreement, 0, 1)';

%% displaying macronutrient agreement
disp('--- RESULTS ---')
for mm = 1:length(str)
    disp(['Macronutrient agreement for ' str{mm} ': ' num2str(agreement(mm)*100) '% | SD: ' num2str(std_dev_performance(mm))])
end

disp(' ')
disp('Press any key to continue')
pause()

agreement = table(str,agreement,std_dev_performance,'VariableNames',{'Macronutrient type','Mean agreement','Standard Deviation'});

end
