function [idxk,bestk,meanS,sdS] = kmeans_bestk(params)
%KMEANS_BESTK performs k-means clustering of params (data matrix) and finds
%the best k as the one maximizing the mean silhouette, for k = 2:6
% OUTPUT:
% idxk: cluster assignments obtained with bestk
% bestk: best number of clusters
% meanS, sdS: mean and SD of the silhouette values for k = 2:6

rng(1)

%% CLUSTERING
meanS = zeros(1,5); sdS = meanS;
for kk = 2:6  % for each k
    idxk = kmeans(params,kk,"Replicates",5,"Distance","sqeuclidean",'Display','off',MaxIter=1e4);

    S = silhouette(params,idxk,"sqEuclidean");
    meanS(kk-1) = mean(S);
    sdS(kk-1) = std(S,[],1);
end

% best number of clusters = max of the mean silhouette
[~,tmp] = max(meanS);
bestk = tmp+1;

% re-computing k-means clustering with best k
idxk = kmeans(params, bestk, "Replicates",5,"Distance","sqeuclidean",MaxIter=1e4);

end
