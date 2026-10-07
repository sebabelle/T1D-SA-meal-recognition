function DB = DB_index(X, labels)

%DB_INDEX computes the Davies-Bouldin index
% ATTENTION: the global mean is assumed to be 0, hence parset must be
% z-scored (using the same observations that are clustered).
% INPUT:
% parset: matrix of z-scored parameters (#observations x #parameters)
% idxk: cluster assignment

clusters = unique(labels);
k = length(clusters);
[n,d] = size(X);

centroids = zeros(k,d);
S = zeros(k,1);

% Compute centroids
for i = 1:k
    Ci = X(labels==clusters(i),:);
    centroids(i,:) = mean(Ci);
end

% Compute intra-cluster scatter
for i = 1:k
    Ci = X(labels==clusters(i),:);
    S(i) = mean(sqrt(sum((Ci - centroids(i,:)).^2,2)));
end

% Compute centroid distances
M = pdist2(centroids,centroids);

R = zeros(k,k);

for i = 1:k
    for j = 1:k
        if i ~= j
            R(i,j) = (S(i) + S(j)) / M(i,j);
        end
    end
end

D = max(R,[],2);

DB = mean(D);

end