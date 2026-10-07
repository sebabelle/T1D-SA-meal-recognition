function [CH_index] = CH_index(parset,idxk)
%CH_INDEX computes the Calinski-Harabasz index
% ATTENTION: the global mean is assumed to be 0, hence parset must be
% z-scored (using the same observations that are clustered).
% INPUT:
% parset: matrix of z-scored parameters (#observations x #parameters)
% idxk: cluster assignment
arguments (Input)
    parset
    idxk
end

arguments (Output)
    CH_index
end

BCSS = 0; % between-cluster sum of squares
WCSS = 0; % within-cluster sum of squares
tot_obs = size(parset,1);
n_clust = length(unique(idxk));

for kk = 1:n_clust
    cluster_mean = mean(parset(idxk == kk, :), 1);
    BCSS = BCSS + sum(cluster_mean.^2) * sum(idxk == kk);
    WCSS = WCSS + sum(sum(pdist2(parset(idxk == kk, :), cluster_mean).^2));
end
CH_index = (BCSS/(n_clust-1))/(WCSS/(tot_obs-n_clust)); 

end
