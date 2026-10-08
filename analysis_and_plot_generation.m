%% GMM CLUSTERING: TRAINING, TEST AND PLOT GENERATION
% Fits the GMM (features and number of clusters selected with
% hyperparameter_tuning.m) on the training subjects, assigns the test
% subjects through the posterior probability, computes the clustering
% metrics and generates the figures.
% Requires: Statistics and Machine Learning Toolbox
% Input (working directory): PARAMETERS_notNormalized.csv, Dataset/ folder, functions/ folder
% Output: assignments.mat, posterior.mat, gmm_model.mat
clear all, close all, clc

addpath("functions")

paramsTable_notNormalized = readtable("PARAMETERS_notNormalized.csv");
paramsTable_notNormalized = paramsTable_notNormalized(:,[6 8 10:34]);

fullTABLE_notNormalized = readtable("PARAMETERS_notNormalized.csv");

% macronutrient variables of interest
str = {'Fat_amt';'Prot_amt'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       CREATE TRAINING TEST SET                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subject-level hold-out (70% training / 30% test)
rng(1)
ID = table2array(fullTABLE_notNormalized(:,"ID"));
uniqueID = unique(ID);
cv = cvpartition(length(uniqueID),'HoldOut',0.3);
trainID = uniqueID(training(cv));
testID = uniqueID(test(cv));
idxTraining = ismember(ID, trainID);
idxTest  = ismember(ID, testID);

paramsTableTest_notNormalized = paramsTable_notNormalized(idxTest,:);
fullTABLETest_notNormalized = fullTABLE_notNormalized(idxTest,:);

paramsTable_notNormalized = paramsTable_notNormalized(idxTraining,:);
paramsTable = normalize(paramsTable_notNormalized,"zscore");
fullTABLE_notNormalized = fullTABLE_notNormalized(idxTraining,:);
fullTABLE = normalize(fullTABLE_notNormalized,"zscore");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           GAUSSIAN MIXTURES WITH FEW PARAMETERS                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp(' --------- Clustering with set of params ---------')
disp(' ')

% selecting features to perform clustering with
parNames = {'RaAUC120min','TGR50','d'};

disp('The set of params used are: ')
parNames

% Extract parameters indeces
allNames = paramsTable_notNormalized.Properties.VariableNames;
[tf, parset_indeces] = ismember(parNames, allNames);

% extracting parameters in a matrix
params_notNormalized = table2array(paramsTable_notNormalized);
[params, mu_zscore, sigma_zscore] = zscore(params_notNormalized);
parset = params(:,parset_indeces); % RaAUC120min, TGR50%, d
parset_notNormalized = params_notNormalized(:,parset_indeces);

rng(1)

%% finding best number of clusters
[n_clust, meanS, sdS] = gmm_bestk(parset,6);

rng(1)
% n_clust = 2; % uncomment to override the gmm_bestk decision

%% Performing unsupervised CV-accuracy
parset_prediction = cv_gm_with_SD(fullTABLE_notNormalized,str,n_clust,parset_notNormalized);
% mean CV accuracy of fat and protein
meanCV = mean(table2array(parset_prediction([1 2],2)));
%%

% performing clustering
gm_dist = fitgmdist(parset,n_clust, ...
    'RegularizationValue',0.01,'CovarianceType','diagonal','SharedCovariance',true);
idxk_train = cluster(gm_dist,parset); % identifying clusters

%% finding most important features according to Mahalanobis
mahalanobisDistance(gm_dist,paramsTable,parset_indeces);
%%

%% statistical test for features,uncomment to show
% for feature = 1:length(parset_indeces)
%     ranksum(parset(idxk_train == 1,feature),parset(idxk_train == 2,feature))
% end
%%

%% Computing entropy
p = gm_dist.ComponentProportion;
entropy = -sum(p.*log2(p))/log2(n_clust);
%%

%% Computing entropy-weighted CV accuracy
wCV = entropy*meanCV;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  CLUSTERING VALIDITY INDICES                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DB = DB_index(parset,idxk_train);
CH = CH_index(parset,idxk_train);

%% calculating numerosity of clusters
for ii = 1:n_clust % for each cluster
    n_observations(ii,:) = sum(idxk_train==ii);
end
%%

%% calculating the median
medianValues = clustersMedian(fullTABLE,str,n_clust,idxk_train);

% identifying labels: FA = lowest median Fat, SA = highest median Fat
[~,idxMedians] = sort(medianValues(1,:),2);
idxFA = idxMedians(1); idxSA = idxMedians(2);
%%

%% Calculating posterior probability
post = posterior(gm_dist,parset);
%%

%% SOFT-CLUSTERING
% observations with a posterior probability within the threshold can be
% assigned to both clusters
threshold = [0.3 0.7];
idxBoth = find(post(:,1)>=threshold(1) & post(:,1)<=threshold(2)); 
numInBoth = length(idxBoth);

disp(['There are ' num2str(numInBoth) ' observations that can be in either clusters'])
disp(' ')

%% STATISTICAL TESTS
macronutrientsMatrix = table2array(fullTABLE(:,str));
MacronutrientsSumRankTest(macronutrientsMatrix,idxk_train)

%% SAVING THE TRAINED MODEL
gmm_model.gm_dist   = gm_dist;
gmm_model.featNames = parNames;
gmm_model.mu        = mu_zscore(parset_indeces);
gmm_model.sigma     = sigma_zscore(parset_indeces);
gmm_model.idxFA      = idxFA;
gmm_model.idxSA      = idxSA;
gmm_model.threshold = threshold;
save gmm_model.mat gmm_model -mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         COLORS                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB default color order (first 2 colors are used by boxchart)
defaultColors = get(groot,'defaultAxesColorOrder');
orange = defaultColors(2,:);   % cluster SA color
blue   = defaultColors(1,:);   % cluster FA color

% colormap with a smooth gradient between the two colors
nShades = 256;
cmap = flipud([linspace(blue(1), orange(1), nShades)', ...
    linspace(blue(2), orange(2), nShades)', ...
    linspace(blue(3), orange(3), nShades)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           (a) SCATTER: CLUSTERING RESULTS (TRAINING)                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
tiledlayout(3,2,"TileSpacing","compact","Padding","compact")

nexttile
scatter3(parset(:,3),parset(:,2),parset(:,1),30,post(:,idxFA),'filled')
hold on

zlabel('AURa(120)','Interpreter','tex')
ylabel('TGR50','Interpreter','tex')
xlabel('d','Interpreter','tex')

view([292.00 18.00])
xlim([-4 4]); ylim([-2 3]); zlim([-2 3])

colormap(cmap)
c2 = colorbar; 
alpha(0.5)     % transparency of dots
title('(a) Clustering results','FontSize',14)
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  (b) BOXCHART: LABELLING RESULTS (TRAINING)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

% Fat and Protein (not normalized), stacked, with cluster as grouping variable
data = table2array(fullTABLE_notNormalized(:,{'Fat_amt','Prot_amt'}));
xPlot = [ones(length(data),1); 2*ones(length(data),1)]; % 1 = Fat, 2 = Protein
allData = [data(:,1); data(:,2)];
allGroups = [idxk_train; idxk_train];

b = boxchart(xPlot, allData, 'GroupByColor', allGroups,'LineWidth',1.4);
b(1).BoxWidth = 0.8; b(2).BoxWidth = 0.8;
b(idxSA).BoxFaceColor = orange; b(idxFA).BoxFaceColor = blue;
b(idxSA).BoxEdgeColor = orange; b(idxFA).BoxEdgeColor = blue;
b(idxSA).MarkerColor = orange; b(idxFA).MarkerColor = blue;

xlim([0.5, 2 + 0.5])

% significance asterisks (positions set manually)
xAsterisk = [1, 2];
yAsterisk = [140 140];
hold on
text(xAsterisk, yAsterisk, '*', 'HorizontalAlignment','center', 'FontSize',18)

xticks(1:2)
xticklabels({'Fat','Protein'})
ylim([0 180])
ax = gca;
ax.FontSize = 12;
axis("square")

ylabel('amount [g]')
legend(b([idxFA idxSA]),{'FA','SA'})
title('(b) Labelling results','FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                (c) BOXCHART: FEATURE DISTRIBUTION (TRAINING)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(data);

% long-form table for boxchart
featNames = {'AURa(120min)','TGR50','d'};
Nfeat = length(featNames);
G = repelem((1:Nfeat)',N,1);          % feature indicator repeated for each observation
Y = reshape(parset,[],1);             % stacked values (Nfeat*N x 1)
Group = repmat(idxk_train,Nfeat,1);   % grouping variable repeated for each feature
Feature = categorical(G,1:Nfeat,featNames);

T = table(Y,Feature,Group);

nexttile(3,[1 2])
hold on;

% cluster colors: row g = color of cluster g
colors = zeros(2,3);
colors(idxFA,:) = blue; colors(idxSA,:) = orange;

% boxcharts of each group with horizontal offsets
groups = unique(Group);
offset = 0.35; % horizontal offset for separating group boxes
hBox = gobjects(numel(groups),1);
for i = 1:numel(groups)
    g = groups(i);
    sel = T.Group==g;
    xpos = double(T.Feature(sel));
    xpos = xpos + (i-1.5)*offset;

    hBox(g) = boxchart(xpos, T.Y(sel), 'BoxWidth', 0.25,'LineWidth',1.4);
    hBox(g).BoxFaceColor = colors(g,:);
    hBox(g).MarkerColor  = colors(g,:);
end
xline(1.5,'Color',[0 0 0 0.5],'LineWidth',0.75)
xline(2.5,'Color',[0 0 0 0.5],'LineWidth',0.75)
xline(3.5,'Color',[0 0 0 0.5],'LineWidth',0.75)
yline(5,'Color',[0 0 0 0.5],'LineWidth',0.75)

set(gca, 'XTick', 1:Nfeat), 
ax = gca;
ax.FontSize = 12;
ax.XTickLabels= {'AURa(120)', ...
    'TGR50', ...
    'd'};
ax.TickLabelInterpreter = 'tex';

ylim([-4 5])
xlim([0.5 3.5]);
ylabel('normalized value');
legend(hBox([idxFA idxSA]),{'FA','SA'}, 'Location', 'south');
title('(c) Feature distributions','FontSize',14)

%% Not normalized features of interest (mean and SD of each cluster)
parset_table_notNormalized = paramsTable_notNormalized(:,{'RaAUC120min','TGR50','d'});

pars_SA = parset_table_notNormalized{idxk_train == idxSA, :};
pars_FA = parset_table_notNormalized{idxk_train == idxFA, :};

stats_pars_SA = [mean(pars_SA); std(pars_SA, [], 1)];
stats_pars_FA = [mean(pars_FA); std(pars_FA, [], 1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  CGM TRACKS: TRAINING SET                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SA_meals_idx = find(post(:,idxSA)>0.5);
FA_meals_idx = find(post(:,idxSA)<0.5);

[t_CGM_grid,~,m_SA,l_p_SA,u_p_SA] = extractingCGMtracks(SA_meals_idx,idxTraining,'SA');
[~,~,m_FA,l_p_FA,u_p_FA] = extractingCGMtracks(FA_meals_idx,idxTraining,'FA');

plotMedianCGM(t_CGM_grid,m_SA,l_p_SA,u_p_SA,m_FA,l_p_FA,u_p_FA,orange,blue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TEST SET                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('------- Moving to test set ---------')
pause()

%% NORMALIZING TEST SET (with the training mean and SD)
% parameters
paramsTest_notNormalized = table2array(paramsTableTest_notNormalized);
paramsTest = (paramsTest_notNormalized - mu_zscore)./sigma_zscore;
parset_test = paramsTest(:,parset_indeces);

% macronutrients
macronutrients_notNormalized = table2array(fullTABLE_notNormalized(:,str));
[~,mu_zscore_macros,sigma_zscore_macros] = zscore(macronutrients_notNormalized);

macronutrientsTest_notNormalized = table2array(fullTABLETest_notNormalized(:,str));
macronutrientsMatrixTest = (macronutrientsTest_notNormalized - mu_zscore_macros)./sigma_zscore_macros;

%% LOG-LIKELIHOOD
avgLogLtest = mean(log(pdf(gm_dist, parset_test)));
avgLogLtrain = mean(log(pdf(gm_dist, parset)));

%% ASSIGNMENT AND STATISTICAL TEST
testPost = posterior(gm_dist,parset_test);
[~,idxk_test] = max(testPost,[],2); % predicted label

MacronutrientsSumRankTest(macronutrientsMatrixTest,idxk_test)

%% NUMBER OF OBSERVATIONS
for ii = 1:n_clust % for each cluster
    n_observations_test(ii,:) = sum(idxk_test==ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               (d) SCATTER: ESTIMATED POSTERIOR (TRAINING + TEST)        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

nexttile
scatter3(parset(:,3),parset(:,2),parset(:,1),30,post(:,idxFA),'LineWidth',1.2)
alpha(0.3)
hold on
scatter3(parset_test(:,3),parset_test(:,2),parset_test(:,1),30,testPost(:,idxFA),'filled')
hold on
view([292.00 18.00])
xlim([-4 4]); ylim([-2 3]); zlim([-2 3])

zlabel('AURa(120)','Interpreter','tex')
ylabel('TGR50','Interpreter','tex')
xlabel('d','Interpreter','tex')

colormap(cmap)
c2 = colorbar; 
title('(d) Estimated posterior','FontSize',14)
legend('Training','Test','Location','northwest')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  (e) BOXCHART: ASSIGNED LABELS (TEST)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

% Fat and Protein (not normalized), stacked, with cluster as grouping variable
data = table2array(fullTABLETest_notNormalized(:,{'Fat_amt','Prot_amt'}));
xPlot = [ones(length(data),1); 2*ones(length(data),1)]; % 1 = Fat, 2 = Protein
allData = [data(:,1); data(:,2)];
allGroups = [idxk_test; idxk_test];

b = boxchart(xPlot, allData, 'GroupByColor', allGroups,'LineWidth',1.4);
b(1).BoxWidth = 0.8; b(2).BoxWidth = 0.8;
b(idxSA).BoxFaceColor = orange; b(idxFA).BoxFaceColor = blue;
b(idxSA).BoxEdgeColor = orange; b(idxFA).BoxEdgeColor = blue;
b(idxSA).MarkerColor = orange; b(idxFA).MarkerColor = blue;

xlim([0.5, 2 + 0.5])

% significance asterisks (positions set manually)
xAsterisk = [1, 2];
yAsterisk = [140 140];
hold on
text(xAsterisk, yAsterisk, '*', 'HorizontalAlignment','center', 'FontSize',18)

xticks(1:2)
xticklabels({'Fat','Protein'})
ylim([0 180])
ax = gca;
ax.FontSize = 12;
axis("square")

ylabel('amount [g]')
legend(b([idxFA idxSA]),{'FA','SA'})
title('(e) Assigned labels','FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  CGM PROFILES: TEST SET                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SA_meals_idx = find(testPost(:,idxSA)>0.5);
FA_meals_idx = find(testPost(:,idxSA)<0.5);

[t_CGM_grid,~,m_SA,l_p_SA,u_p_SA] = extractingCGMtracks(SA_meals_idx,idxTest,'SA');
[~,~,m_FA,l_p_FA,u_p_FA] = extractingCGMtracks(FA_meals_idx,idxTest,'FA');

plotMedianCGM(t_CGM_grid,m_SA,l_p_SA,u_p_SA,m_FA,l_p_FA,u_p_FA,orange,blue)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  CGM PROFILES: TRAINING + TEST                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAINING => clustering results
% TEST => assignment results
N_meals = length(idxTraining);
idxk = zeros(N_meals,1);
idxk(idxTraining) = idxk_train;
idxk(idxTest) = idxk_test;

SA_meals_idx = find(idxk == idxSA);
FA_meals_idx = find(idxk == idxFA); 

[t_CGM_grid,~,m_SA,l_p_SA,u_p_SA] = extractingCGMtracks(SA_meals_idx,ones(N_meals,1),'SA');
[~,~,m_FA,l_p_FA,u_p_FA] = extractingCGMtracks(FA_meals_idx,ones(N_meals,1),'FA');

plotMedianCGM(t_CGM_grid,m_SA,l_p_SA,u_p_SA,m_FA,l_p_FA,u_p_FA,orange,blue)


%% ============= CREATE ASSIGNMENTS ============= %%
% 1 = meal assigned to cluster SA, 0 = meal assigned to cluster FA
assignments = zeros(N_meals,1);

assignments(idxTraining) = idxk_train == idxSA;
assignments(idxTest) = idxk_test == idxSA;

save assignments.mat assignments -mat

%% =============== CREATE POSTERIOR ============== %%
trainPost = post;

post = zeros(N_meals,2);

post(idxTraining,:) = trainPost;
post(idxTest,:) = testPost;

save posterior.mat post -mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          LOCAL FUNCTIONS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMedianCGM(t_grid,m_SA,l_p_SA,u_p_SA,m_FA,l_p_FA,u_p_FA,orange,blue)
    % Median CGM profile (and 25th-75th percentile band) of cluster SA and FA
    figure
    % SA
    plot(t_grid,m_SA,'Color',orange,'LineWidth',2)
    hold on
    fill([t_grid; flipud(t_grid)], [l_p_SA; flipud(u_p_SA)], ...
        orange,'FaceAlpha',0.5,'EdgeColor','none')
    % FA
    plot(t_grid,m_FA,'Color',blue,'LineWidth',2)
    fill([t_grid; flipud(t_grid)], [l_p_FA; flipud(u_p_FA)], ...
        blue,'FaceAlpha',0.5,'EdgeColor','none')
    xline(0,'LineWidth',1,'Color','k','LineStyle','-.') % meal time
    xlabel('time [min]'), ylabel('CGM [mg/dl]')
    grid on
    yline(180,'LineWidth',1,'LineStyle','--')
    yline(70,'LineWidth',1,'LineStyle','--')
    legend('','SA','','FA','','Location','southeast')
    xlim([-60 240])
    ylim([60 240])
end
