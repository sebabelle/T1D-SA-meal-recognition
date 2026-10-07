function [idxk,bestk,meanS,sdS] = hclust_bestk(params)
%HCLUST_BESTK performs hierarchical clustering (Ward linkage) of params
%(data matrix) and finds the best k as the one maximizing the mean
%silhouette, for k = 2:6
% OUTPUT:
% idxk: cluster assignments obtained with bestk
% bestk: best number of clusters
% meanS, sdS: mean and SD of the silhouette values for k = 2:6

rng(1)

%% creating the tree with Ward's linkage method
h_tree = linkage(params,"ward");

%% CLUSTERING
for kk = 2:6  % for each k
    idxk = cluster(h_tree,"maxclust",kk);

    S = silhouette(params,idxk,"sqEuclidean");
    meanS(kk-1) = mean(S);
    sdS(kk-1) = std(S,[],1);
end

% best number of clusters = max of the mean silhouette
[~,tmp] = max(meanS);
bestk = tmp+1;

idxk = cluster(h_tree,"maxclust",bestk);

end
