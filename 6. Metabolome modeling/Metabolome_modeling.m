% PARAFAC functionality
home = ".";
cd(home)
addpath("..\Matlab scripts\Scripts\"); % own scripts
addpath("..\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
% Load input data
numTimepoints_metab = 5;

metabolomics_cnt_scl = readmatrix("../3. Metabolome pre-processing/Metabolomics.csv");
metabolomics_cnt_scl = reshape(metabolomics_cnt_scl, size(metabolomics_cnt_scl, 1), size(metabolomics_cnt_scl,2)/numTimepoints_metab, numTimepoints_metab);

%%
% Load metadata
metabolomics_id_meta = readmatrix("../3. Metabolome pre-processing/Metabolomics_id_meta.csv", OutputType="string");
metabolomics_feature_meta = readmatrix("../3. Metabolome pre-processing/Metabolomics_feature_meta.csv", OutputType="string");

%%
% Import red fluorescence data
rf_data = readmatrix("../0. Raw data input/RFdata.csv", OutputType="string");
rf_data = rf_data(:, [1 6]);     % keep subject + RF group information

% Fix incorrect subject names
rf_data(rf_data(:,1) == "VSTPHZ", 1) = "VSTPH2";
rf_data(rf_data(:,1) == "D2VZH0", 1) = "DZVZH0";
rf_data(rf_data(:,1) == "DLODNN", 1) = "DLODDN";
rf_data(rf_data(:,1) == "O3VQFX", 1) = "O3VQFQ";
rf_data(rf_data(:,1) == "F80LGT", 1) = "F80LGF";
rf_data(rf_data(:,1) == "26QQR0", 1) = "26QQrO";

rf_data = unique(rf_data, "rows");
rf_data = rf_data(ismember(rf_data(:,1), metabolomics_id_meta),:);

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

%%
% maxComponents=5;
% numReps=5;
% numSplits=10;
% path = "./Figures/Metabolomics/";
% 
% numComponentsOverviewPlot(metabolomics_cnt_scl,maxComponents,numReps,numSplits,"Metabolomics",path+"components_metabolomics.jpg");

%%
% Do a reporter Parafac + plot it to check what the best model is.
maxComponents=3;
path_start = "./20230605_run/Figures/";
days = [0 2 5 9 14];
numReps=100;

[metabolomicsModels, metabolomicsCons, metabolomicsVarExps, metabolomicsBoots, metabolomicsBootsVarExp, metabolomicsTuckers] = quickReport(metabolomics_cnt_scl, maxComponents, numReps, rf_data, metabolomics_feature_meta, days, 1, "Metabolomics bootstrapped", path_start+"Metabolomics");

%%
% Dump data so far for later inspection
path_start = "./20230605_run/Dump/";
dump(metabolomicsModels, metabolomicsCons, metabolomicsVarExps, metabolomicsBoots, metabolomicsBootsVarExp, metabolomicsTuckers, path_start, "metabolomics");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.

numFactors = 2;
%metabolomics_choice = find(metabolomicsVarExps{numFactors}==max(metabolomicsVarExps{numFactors}));
%metabolomics_choice = find(metabolomicsCons{numFactors}==max(metabolomicsCons{numFactors}));
metabolomics_choice = 11; % override because varexps cant find the right one

Metabolomics_model = pickModel(metabolomicsModels{1,numFactors}, metabolomicsModels{2,numFactors}, metabolomicsModels{3,numFactors}, metabolomics_choice);

%%
% Save the models
model_path = "./20230605_run/PARAFAC models/";
savePARAFAC(metabolomics_cnt_scl, Metabolomics_model, rf_data, metabolomics_feature_meta, model_path + "Metabolomics");

%%
% Plot PARAFAC models
days = [0 2 5 9 14];
timepoints = 1:5;
path_start = "./20230605_run/Figures/";

plotPARAFAC4(metabolomics_cnt_scl, Metabolomics_model, metabolomicsVarExps{numFactors}, metabolomics_choice, rf_data, metabolomics_feature_meta, days, timepoints, 1, "PARAFAC metabolomics", path_start + "PARAFAC_metabolomics.jpg");

%%
% Make biplots
path_start = "./20230605_run/Figures/";
biplotPARAFAC(metabolomics_cnt_scl, Metabolomics_model, rf_data, metabolomics_feature_meta, 1, "PARAFAC metabolomics", path_start+"Metabolomics_biplot.jpg");