function [mahalD] = mahalanobisDistance(gm,paramsTable,parset_indeces)
%MAHALANOBISDISTANCE computes, for each parameter, the distance between the
%means of two gaussian components, scaled by the variance
% ATTENTION: designed for a GMM with 2 components
arguments (Input)
   gm % gaussian mixture model object
   paramsTable % table with parameters
   parset_indeces % indices of the parameters used in the GMM
end

arguments (Output)
    mahalD
end

% parameters names
parNames = paramsTable.Properties.VariableNames(parset_indeces);

Ncomponents = gm.NumComponents;

mu = gm.mu;
Sigma = gm.Sigma;

for ii = 1:Ncomponents
    for jj = ii+1:Ncomponents
            diff = mu(ii,:) - mu(jj,:);

            % Mahalanobis distance
            d = sqrt((diff.^2) ./ Sigma);
    end
end

mahalD = table(d','RowNames',parNames','VariableNames',{'Mahalanobis Distance'});

% displaying
mahalD

end
