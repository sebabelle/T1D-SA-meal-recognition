function plotMixtures(gmm_model,Xz)
%% Plot: fitted mixture components and new meals
% Top row: marginal density of each feature (components weighted by their
% mixing proportion, dashed = whole mixture). Bottom row: pairs of features
% with the 1 SD and 2 SD contours of each component. Axes are in z-score
% units of the training set; dots are the new meals colored by posterior FA.

arguments (Input)
    gmm_model % Gaussian Mixture model
    Xz % z-scored observations
end




gm = gmm_model.gm_dist;
K = gm.NumComponents;
nF = size(gm.mu,2);
idxFA = gmm_model.idxFA;
idxSA = gmm_model.idxSA;
defaultColors = get(groot,'defaultAxesColorOrder');
orange = defaultColors(2,:);   % cluster SA color
blue   = defaultColors(1,:);   % cluster FA color
nShades = 256;
cmap = flipud([linspace(blue(1), orange(1), nShades)', ...
    linspace(blue(2), orange(2), nShades)', ...
    linspace(blue(3), orange(3), nShades)']);


compColor = zeros(K,3);
compColor(idxFA,:) = blue; compColor(idxSA,:) = orange;

% posterior
post = posterior(gmm_model.gm_dist,Xz);
postFA = post(:,idxFA);


% standard deviation of each component (diagonal covariance)
sd = zeros(K,nF);
for k = 1:K
    sd(k,:) = sqrt(gm.Sigma); % shared covariance
end

xg = linspace(min([-4; Xz(:)]),max([4; Xz(:)]),400);
pairs = nchoosek(1:nF,2);
nCols = max(nF,size(pairs,1));

figure
tiledlayout(1,nCols,"TileSpacing","compact","Padding","compact")

% pairs of features with 1 SD and 2 SD limits
theta = linspace(0,2*pi,200);
lineStyles = {'-','--'};
for pp = 1:size(pairs,1)
    a = pairs(pp,1); c = pairs(pp,2);
    ax = nexttile(pp);
    hold on
    scatter(Xz(:,a),Xz(:,c),25,postFA,'filled','MarkerFaceAlpha',0.6)
    for k = 1:K
        for r = 1:2
            plot(gm.mu(k,a) + r*sd(k,a)*cos(theta), gm.mu(k,c) + r*sd(k,c)*sin(theta), ...
                'Color',compColor(k,:),'LineWidth',1.5,'LineStyle',lineStyles{r})
        end
        plot(gm.mu(k,a),gm.mu(k,c),'x','Color',compColor(k,:),'LineWidth',2,'MarkerSize',10)
    end
    colormap(ax,cmap)
    ax.CLim = [0 1];
    xlabel(gmm_model.featNames{a},'Interpreter','none')
    ylabel(gmm_model.featNames{c},'Interpreter','none')
    grid on
    if pp == size(pairs,1)
        cb = colorbar(ax); cb.Label.String = 'posterior FA';
    end
end

end