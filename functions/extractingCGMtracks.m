function [t_grid, mean_curve, median_curve, low_prc, up_prc] = extractingCGMtracks(idx_meals_to_extract, idxcvpartition, meal_label)
%EXTRACTINGCGMTRACKS extracts the time and CGM vectors of each meal in
%idx_meals_to_extract, plots each CGM track together with the median curve
%and the 25th-75th percentile band.
%% ATTENTION: the Dataset folder (one .mat file per meal) must be in the working directory

arguments (Input)
    idx_meals_to_extract % indices of the meals to be extracted (e.g. the meals of a cluster)
    idxcvpartition       % logical/index vector selecting the meals of the subset analysed (e.g. training set)
end

arguments (Input, Repeating)
    meal_label % OPTIONAL: label of the analysed meals (char), used in the plot title
end

arguments (Output)
    t_grid       % bin centers [min]
    mean_curve   % mean curve (bins of 10 min width)
    median_curve % median curve (bins of 10 min width)
    low_prc      % lower percentile (25%)
    up_prc       % upper percentile (75%)
end

%% Checking existence of the Dataset folder
if ~isfolder('Dataset')
    error('extractingCGMtracks:DatasetNotFound','Dataset folder not found in the working directory')
end
addpath('Dataset')
dataset_info = dir(fullfile('Dataset','*.mat'));

%% Loading CGM tracks
n_meals_to_extract = numel(idx_meals_to_extract);

figure
all_t_CGM = []; all_CGM = [];
idxcvpartition = find(idxcvpartition);
for jj = 1:n_meals_to_extract % for each meal to extract
    data = load(dataset_info(idxcvpartition(idx_meals_to_extract(jj))).name);

    t_CGM = data.data.t_CGM;
    all_t_CGM = [all_t_CGM; t_CGM];

    CGM = data.data.CGM;
    all_CGM = [all_CGM; CGM];

    % spaghetti plot with white underlay for readability
    plot(t_CGM,CGM,'Color',[1 1 1],'LineWidth',1)
    hold on
    plot(t_CGM,CGM,'Color',[0.8 0.8 0.8])
end

%% Binning
binWidth = 10/60;  % hours
binEdges = min(all_t_CGM/60):binWidth:max(all_t_CGM/60);
binCenters = (binEdges(1:end-1) + binWidth/2)';

% adding the end points (-60 and 240 min) so that the curves span the whole plot
binCenters = [-1; binCenters; 4]*60; % [min]

[median_curve,low_prc,up_prc,mean_curve,~] = ...
    prctiles_curve(all_CGM,all_t_CGM/60,binCenters/60,25,75);

t_grid = binCenters;

%% Plotting 25th-75th percentiles and median
fill([binCenters; flipud(binCenters)],[low_prc; flipud(up_prc)],...
    'k','FaceAlpha',0.3,'EdgeColor','none')
hold on
plot(t_grid,median_curve,'LineStyle','-','LineWidth',2,'Color','k')

xlabel('time [min]')
ylabel('CGM [mg/dl]')
xlim([-60 240])
grid on

if ~isempty(meal_label) && ischar(meal_label{1})
    title(['CGM curves of meals belonging to cluster ' meal_label{1}])
else
    title('CGM curves')
end

end
