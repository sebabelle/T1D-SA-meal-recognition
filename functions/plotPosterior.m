function [outputArg1,outputArg2] = plotPosterior(gmm_model,Xz,label)
%% Plot: estimated posterior in the 3D feature space

arguments (Input)
    gmm_model % Gaussian Mixture model
    Xz % z-scored observations
    label % assigned labels
end

post = posterior(gmm_model.gm_dist,Xz);

idxFA = gmm_model.idxFA; idxSA = gmm_model.idxSA;
postFA = post(:,idxFA);

defaultColors = get(groot,'defaultAxesColorOrder');
orange = defaultColors(2,:);   % cluster SA color
blue   = defaultColors(1,:);   % cluster FA color
nShades = 256;
cmap = flipud([linspace(blue(1), orange(1), nShades)', ...
    linspace(blue(2), orange(2), nShades)', ...
    linspace(blue(3), orange(3), nShades)']);

figure
scatter3(Xz(:,3),Xz(:,2),Xz(:,1),50,postFA,'filled')
text(Xz(:,3),Xz(:,2),Xz(:,1),label,"HorizontalAlignment","center",'VerticalAlignment','bottom')
xlabel(gmm_model.featNames{3},'Interpreter','none')
ylabel(gmm_model.featNames{2},'Interpreter','none')
zlabel(gmm_model.featNames{1},'Interpreter','none')
view([292.00 18.00])
colormap(cmap), colorbar
grid on
title('Estimated posterior')


end