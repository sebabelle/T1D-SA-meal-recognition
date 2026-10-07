function [best_k, meanS, sdS] = gmm_bestk(parset,max_clust)
%GMM_BESTK finds the best number of gaussian components (clusters) using
%the silhouette method. The silhouette is computed on the hard-assignment
%based on the maximum posterior.
arguments (Input)
    parset % matrix of parameters #n x pars 
    max_clust % maximum number of gaussian components to test
end

arguments (Output)
    best_k % best number of gaussian components 
    meanS  % mean silhouette values for each number of clusters tested (2:max_clust)
    sdS    % standard deviation of the silhouette for each number of clusters tested
end

rng(1)

%% Defining best number of clusters
meanS = zeros(max_clust-1,1); sdS = meanS;
for n_clust = 2:max_clust
    curr_gm_dist = fitgmdist(parset,n_clust, ...
    'RegularizationValue',0.01,'CovarianceType','diagonal','SharedCovariance',true);

    % hard-assignment based on the component with the maximum posterior
    curr_post = posterior(curr_gm_dist,parset);
    [~,curr_idx] = max(curr_post,[],2);

    meanS(n_clust-1,:) = mean(silhouette(parset,curr_idx));
    sdS(n_clust-1,:) = std(silhouette(parset,curr_idx),[],1);
end

% maximum of the mean silhouette
[~,best_k] = max([0; meanS]);

end
