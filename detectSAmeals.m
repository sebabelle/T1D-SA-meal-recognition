function results = detectSAmeals(data,options)
%CLASSIFYMEALS assigns meals to cluster FA or SA using the trained GMM
%(gmm_model.mat). New data are normalized with the mean and SD of the
%training set.
% USAGE:
% results = detectSAmeals("meals.csv")
% results = detectSAmeals(T,ModelFile="gmm_model.mat",OutputFile="out.csv", ...
%                         PlotPosterior=true,PlotMixtures=true)
%           where T is a MATLAB table
% INPUT:
% data: path of a CSV file or a table, one row per meal, with the (not
%       normalized) columns RaAUC120min, TGR50, d. An "ID" column is optional.
% OPTIONAL 
% ModelFile: trained model (default: gmm_model.mat in the folder of this function)
% OutputFile: if given, the results are saved as CSV
% PlotPosterior: 3D scatter of the meals colored by posterior FA (default false)
% PlotMixtures: fitted mixture components with the meals on top (default false)
% OUTPUT:
% results: table with Posterior_FA, Posterior_SA, Label (FA/SA, maximum
%          posterior) and Uncertain (posterior SA within the model thresholds).
%          Meals with missing values are not assigned (NaN / empty label).
% Requires: Statistics and Machine Learning Toolbox
arguments
    data
    options.ModelFile (1,1) string = fullfile(pwd,"model/gmm_model.mat")
    options.OutputFile (1,1) string = ""
    options.PlotPosterior (1,1) logical = false
    options.PlotMixtures (1,1) logical = false
end

%% Loading model and data
if ~isfile(options.ModelFile) % if model/gmm_model.mat does not exist
    error('%s not found. Make sure it is in working derectory' , options.ModelFile)
end
load(options.ModelFile,"gmm_model")

% check if input is MATLAB table or .csv
if istable(data)
    newTable = data;
else
    newTable = readtable(data);
end

% check if all features are in table
missingVars = setdiff(gmm_model.featNames, newTable.Properties.VariableNames);
if ~isempty(missingVars)
    error('Columns missing in the input data: %s', strjoin(missingVars,', '))
end

%% Normalizing with the training mean and SD
% z-scoring with mean and SD used in training phase
X = table2array(newTable(:,gmm_model.featNames)); % convert in array
if ~isnumeric(X) % if columns are in not allowed type
    error('The columns %s must be numbers', strjoin(gmm_model.featNames,', '))
end
Xz = (X - gmm_model.mu)./gmm_model.sigma; % z-score element-wise

%% Posterior probability 
% meals with missing values are not assigned
nMeals = size(Xz,1);
valid = ~any(isnan(Xz),2); % check if any 
post = nan(nMeals,gmm_model.gm_dist.NumComponents);
% compute posterior of valid entries
post(valid,:) = posterior(gmm_model.gm_dist,Xz(valid,:));

idxFA = gmm_model.idxFA; idxSA = gmm_model.idxSA;
postFA = post(:,idxFA);
postSA = post(:,idxSA);

%% Label assignment (maximization of posterior) and uncertain meals
label = strings(nMeals,1);
isSA = false(nMeals,1);
% maximizing posterior
[~,idxk_new] = max(post(valid,:),[],2);
isSA(valid) = idxk_new == idxSA;
label(valid) = "FA";
label(isSA) = "SA";
% checking if any meal has uncertain posterior
uncertain = postSA >= gmm_model.threshold(1) & postSA <= gmm_model.threshold(2);

results = table(postFA,postSA,label,uncertain, ...
    'VariableNames',{'Posterior_FA','Posterior_SA','Label','Uncertain'});
if ismember('ID',newTable.Properties.VariableNames) % if there is field ID
    results = [newTable(:,'ID') results]; % save results with field
end
results

% output
fprintf('%d meals assigned to SA, %d to FA, %d with uncertain posterior, %d not assigned (missing values)\n', ...
    sum(label=="SA"),sum(label=="FA"),sum(uncertain),sum(~valid))

if options.OutputFile ~= ""
    writetable(results,options.OutputFile)
end

%% Plots

% create colormap of posterior probability
if (options.PlotPosterior || options.PlotMixtures) && isfolder(fullfile(pwd,'functions'))
    addpath(fullfile(pwd,'functions/'))
    defaultColors = get(groot,'defaultAxesColorOrder');
    orange = defaultColors(2,:);   % cluster SA color
    blue   = defaultColors(1,:);   % cluster FA color
    nShades = 256;
    cmap = flipud([linspace(blue(1), orange(1), nShades)', ...
        linspace(blue(2), orange(2), nShades)', ...
        linspace(blue(3), orange(3), nShades)']);
elseif ~isfolder(fullfile(pwd,'functions'))
    error('Functions folder is not in the working directory.');
end

if options.PlotPosterior
    plotPosterior(gmm_model,Xz,label,cmap)
end
if options.PlotMixtures
    plotMixtures(gmm_model,Xz,cmap)
end

end

