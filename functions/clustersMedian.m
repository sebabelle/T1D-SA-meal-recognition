function M = clustersMedian(fullTABLE,str,bestk,idxk)
% clustersMedian calculates the median value of each macronutrient in each cluster
% INPUT:
% fullTABLE: table containing the (z-scored) dataset
% str: cell array of macronutrient names of interest
% bestk: number of clusters
% idxk: cluster assignment indices
% OUTPUT:
% M: #macronutrients x #clusters

M = zeros(length(str),bestk);

for jj = 1:length(str) % for each macronutrient transformation
    currMacro = table2array(fullTABLE(:,str{jj}));
    for ii = 1:bestk % for each cluster
        M(jj,ii) = median(currMacro(idxk == ii));
    end
end
end
