% PARAFAC functionality
home = ".";
cd(home)
addpath("..\Matlab scripts\Scripts\"); % own scripts
addpath("..\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
% Load input data
numTimepoints_microb = 7;

microb_tongue_cnt_scl = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Microb_tongue_functional.csv");
microb_tongue_cnt_scl = reshape(microb_tongue_cnt_scl, size(microb_tongue_cnt_scl, 1), size(microb_tongue_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_lowling_cnt_scl = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Microb_lowling_functional.csv");
microb_lowling_cnt_scl = reshape(microb_lowling_cnt_scl, size(microb_lowling_cnt_scl, 1), size(microb_lowling_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_lowinter_cnt_scl = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Microb_lowinter_functional.csv");
microb_lowinter_cnt_scl = reshape(microb_lowinter_cnt_scl, size(microb_lowinter_cnt_scl, 1), size(microb_lowinter_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_upling_cnt_scl = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Microb_upling_functional.csv");
microb_upling_cnt_scl = reshape(microb_upling_cnt_scl, size(microb_upling_cnt_scl, 1), size(microb_upling_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_upinter_cnt_scl = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Microb_upinter_functional.csv");
microb_upinter_cnt_scl = reshape(microb_upinter_cnt_scl, size(microb_upinter_cnt_scl, 1), size(microb_upinter_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_saliva_cnt_scl = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Microb_saliva_functional.csv");
microb_saliva_cnt_scl = reshape(microb_saliva_cnt_scl, size(microb_saliva_cnt_scl, 1), size(microb_saliva_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

%%
% Load metadata
microb_tongue_id_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Tongue_id_meta.csv", OutputType="string");
microb_lowling_id_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Lowling_id_meta.csv", OutputType="string");
microb_lowinter_id_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Lowinter_id_meta.csv", OutputType="string");
microb_upling_id_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Upling_id_meta.csv", OutputType="string");
microb_upinter_id_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Upinter_id_meta.csv", OutputType="string");
microb_saliva_id_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Saliva_id_meta.csv", OutputType="string");

microb_tongue_feature_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Tongue_feature_meta.csv", OutputType="string");
microb_lowling_feature_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Lowling_feature_meta.csv", OutputType="string");
microb_lowinter_feature_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Lowinter_feature_meta.csv", OutputType="string");
microb_upling_feature_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Upling_feature_meta.csv", OutputType="string");
microb_upinter_feature_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Upinter_feature_meta.csv", OutputType="string");
microb_saliva_feature_meta = readmatrix("../4. Functional microbiome pre-processing/20230618_run/Saliva_feature_meta.csv", OutputType="string");

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
rf_data = rf_data(ismember(rf_data(:,1), microb_tongue_id_meta),:);

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
% path = "./Figures/Functional analysis/Mapped 100/"; %Pathways/"+pathwayOfInterest+"/";
% 
% numComponentsOverviewPlot(microb_tongue_cnt_scl,maxComponents,numReps,numSplits,"Tongue",path+"components_tongue_func.jpg");
% numComponentsOverviewPlot(microb_low_ling_cnt_scl,maxComponents,numReps,numSplits,"Lower jaw, lingual",path+"components_lowling_func.jpg");
% numComponentsOverviewPlot(microb_low_inter_cnt_scl,maxComponents,numReps,numSplits,"Lower jaw, interproximal",path+"components_lowinter_func.jpg");
% numComponentsOverviewPlot(microb_up_ling_cnt_scl,maxComponents,numReps,numSplits,"Upper jaw, lingual",path+"components_upling_func.jpg");
% numComponentsOverviewPlot(microb_up_inter_cnt_scl,maxComponents,numReps,numSplits,"Upper jaw, interproximal",path+"components_upinter_func.jpg");

%%
% Plot the reporter Parafac output
path_start = "./20230618_run/Figures/";
maxComponents=3;
days = [-14 0 2 5 9 14 21];
numReps=25;
maxIterations=20;

[tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers] = quickReport(microb_tongue_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_tongue_feature_meta, days, 2, "Tongue bootstrapped", path_start+"tongue");
[lowlingModels, lowlingCons, lowlingVarExps, lowlingBoots, lowlingBootVarExps, lowlingTuckers] = quickReport(microb_lowling_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_lowling_feature_meta, days, 2, "Low ling bootstrapped", path_start+"lowling");
[lowinterModels, lowinterCons, lowinterVarExps, lowinterBoots, lowinterBootVarExps, lowinterTuckers] = quickReport(microb_lowinter_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_lowinter_feature_meta, days, 2, "Low inter bootstrapped", path_start+"lowinter");
[uplingModels, uplingCons, uplingVarExps, uplingBoots, uplingBootVarExps, uplingTuckers] = quickReport(microb_upling_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_upling_feature_meta, days, 2, "Up ling bootstrapped", path_start+"upling");
[upinterModels, upinterCons, upinterVarExps, upinterBoots, upinterBootVarExps, upinterTuckers] = quickReport(microb_upinter_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_upinter_feature_meta, days, 2, "Up inter bootstrapped", path_start+"upinter");
[salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers] = quickReport(microb_saliva_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_saliva_feature_meta, days, 2, "Saliva bootstrapped", path_start+"saliva");

%%
% Dump data so far for later inspection
path_start = "./20230618_run/Dump/";

dump(tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers, path_start, "tongue");
dump(lowlingModels, lowlingCons, lowlingVarExps, lowlingBoots, lowlingBootVarExps, lowlingTuckers, path_start, "lowling");
dump(lowinterModels, lowinterCons, lowinterVarExps, lowinterBoots, lowinterBootVarExps, lowinterTuckers, path_start, "lowinter");
dump(uplingModels, uplingCons, uplingVarExps, uplingBoots, uplingBootVarExps, uplingTuckers, path_start, "upling");
dump(upinterModels, upinterCons, upinterVarExps, upinterBoots, upinterBootVarExps, upinterTuckers, path_start, "upinter");
dump(salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers, path_start, "saliva");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.

numFactors_tongue = 2;
numFactors_lowling = 1;
numFactors_lowinter = 2;
numFactors_upling = 2;
numFactors_upinter = 1;
numFactors_saliva = 1;

tongue_choice = find(tongueVarExps{numFactors_tongue}==max(tongueVarExps{numFactors_tongue}));
%tongue_choice = find(tongueCons{numFactors_tongue}==max(tongueCons{numFactors_tongue}));

%lowling_choice = find(lowlingVarExps{numFactors_lowling}==max(lowlingVarExps{numFactors_lowling}));
%lowling_choice = find(lowlingCons{numFactors_lowling}==max(lowlingCons{numFactors_lowling}));
lowling_choice = 3; % due to unique situation neither varexp nor corcondia will find the right model

%lowinter_choice = find(lowinterVarExps{numFactors_lowinter}==max(lowinterVarExps{numFactors_lowinter}));
lowinter_choice = find(lowinterCons{numFactors_lowinter}==max(lowinterCons{numFactors_lowinter}));

upling_choice = find(uplingVarExps{numFactors_upling}==max(uplingVarExps{numFactors_upling}));
%upling_choice = find(uplingCons{numFactors_upling}==max(uplingCons{numFactors_upling}));

upinter_choice = find(upinterVarExps{numFactors_upinter}==max(upinterVarExps{numFactors_upinter}));
%upinter_choice = find(upinterCons{numFactors_upinter}==max(upinterCons{numFactors_upinter}));

saliva_choice = find(salivaVarExps{numFactors_saliva}==max(salivaVarExps{numFactors_saliva}));

Tongue_model = pickModel(tongueModels{1,numFactors_tongue}, tongueModels{2,numFactors_tongue}, tongueModels{3,numFactors_tongue}, tongue_choice);
L_ling_model = pickModel(lowlingModels{1,numFactors_lowling}, lowlingModels{2,numFactors_lowling}, lowlingModels{3,numFactors_lowling}, lowling_choice);
L_inter_model = pickModel(lowinterModels{1,numFactors_lowinter}, lowinterModels{2,numFactors_lowinter}, lowinterModels{3,numFactors_lowinter}, lowinter_choice);
U_ling_model = pickModel(uplingModels{1,numFactors_upling}, uplingModels{2,numFactors_upling}, uplingModels{3,numFactors_upling}, upling_choice);
U_inter_model = pickModel(upinterModels{1,numFactors_upinter}, upinterModels{2,numFactors_upinter}, upinterModels{3,numFactors_upinter}, upinter_choice);
Saliva_model = pickModel(salivaModels{1,numFactors_saliva}, salivaModels{2,numFactors_saliva}, salivaModels{3,numFactors_saliva}, saliva_choice);

%%
% Save the models
model_path = "./20230618_run/PARAFAC models/";

savePARAFAC(microb_tongue_cnt_scl, Tongue_model, rf_data, microb_tongue_feature_meta, model_path +  "Tongue_func");
savePARAFAC(microb_lowinter_cnt_scl, L_inter_model, rf_data, microb_lowinter_feature_meta, model_path +  "Low_inter_func");
savePARAFAC(microb_lowling_cnt_scl, L_ling_model, rf_data, microb_lowling_feature_meta, model_path + "Low_ling_func");
savePARAFAC(microb_upinter_cnt_scl, U_inter_model, rf_data, microb_upinter_feature_meta, model_path + "Up_inter_func");
savePARAFAC(microb_upling_cnt_scl, U_ling_model, rf_data, microb_upling_feature_meta, model_path + "Up_ling_func");
savePARAFAC(microb_saliva_cnt_scl, Saliva_model, rf_data, microb_saliva_feature_meta, model_path + "Saliva_func");

%%
% Plot PARAFAC models
days = [-14 0 2 5 9 14 21];
timepoints = 1:7;
path_start = "./20230618_run/Figures/";

plotPARAFAC3_func(microb_tongue_cnt_scl, Tongue_model, 0, tongueVarExps{numFactors_tongue}, tongue_choice, rf_data, microb_tongue_feature_meta, days, timepoints, "PARAFAC tongue", path_start + "PARAFAC_tongue_func.jpg");
plotPARAFAC3_func(microb_lowinter_cnt_scl, L_inter_model, 0, lowinterVarExps{numFactors_lowinter}, lowinter_choice, rf_data, microb_lowinter_feature_meta, days, timepoints, "PARAFAC lower jaw, interproximal", path_start + "PARAFAC_low_inter_func.jpg");
plotPARAFAC3_func(microb_lowling_cnt_scl, L_ling_model, 0, lowlingVarExps{numFactors_lowling}, lowling_choice, rf_data, microb_lowling_feature_meta, days, timepoints, "PARAFAC lower jaw, lingual", path_start + "PARAFAC_low_ling_func.jpg");
plotPARAFAC3_func(microb_upinter_cnt_scl, U_inter_model, 0, upinterVarExps{numFactors_upinter}, upinter_choice, rf_data, microb_upinter_feature_meta, days, timepoints, "PARAFAC upper jaw, interproximal", path_start + "PARAFAC_up_inter_func.jpg");
plotPARAFAC3_func(microb_upling_cnt_scl, U_ling_model, 0, uplingVarExps{numFactors_upling}, upling_choice, rf_data, microb_upling_feature_meta, days, timepoints, "PARAFAC upper jaw, lingual", path_start + "PARAFAC_up_ling_func.jpg");
plotPARAFAC3_func(microb_saliva_cnt_scl, Saliva_model, 0, salivaVarExps{numFactors_saliva}, saliva_choice, rf_data, microb_saliva_feature_meta, days, timepoints, "PARAFAC saliva", path_start + "PARAFAC_saliva_func.jpg");

%%
% Make biplots
path_start = "./20230618_run/Figures/";

biplotPARAFAC_func(microb_tongue_cnt_scl, Tongue_model, rf_data, microb_tongue_feature_meta, 2, "PARAFAC Tongue", path_start+"Tongue_func_biplot.jpg")
%biplotPARAFAC_func(microb_lowling_cnt_scl, L_ling_model, rf_data, microb_lowling_feature_meta, 2, "PARAFAC Lower jaw, lingual", path_start+"Low_ling_func_biplot.jpg")
biplotPARAFAC_func(microb_lowinter_cnt_scl, L_inter_model, rf_data, microb_lowinter_feature_meta, 2, "PARAFAC Lower jaw, interproximal", path_start+"Low_inter_func_biplot.jpg")
biplotPARAFAC_func(microb_upling_cnt_scl, U_ling_model, rf_data, microb_upling_feature_meta, 2, "PARAFAC Upper jaw, lingual", path_start+"Up_ling_func_biplot.jpg")
%biplotPARAFAC_func(microb_up_inter_cnt_scl, U_inter_model, microb_up_inter_id_meta, microb_up_inter_feature_meta, 2, "PARAFAC Upper jaw, interproximal", path_start+"Up_inter_func_biplot.jpg")

%%
% Check if variation explained is higher than permuted data
% numReps=25;
% 
% permutedPARAFAC(microb_tongue_cnt_scl,2,numReps,"./Permuted PARAFAC models/Tongue_func_2");
% permutedPARAFAC(microb_low_ling_cnt_scl,2,numReps,"./Permuted PARAFAC models/Low_ling_func_2");
% permutedPARAFAC(microb_low_inter_cnt_scl,2,numReps,"./Permuted PARAFAC models/Low_inter_func_2");
% permutedPARAFAC(microb_up_ling_cnt_scl,2,numReps,"./Permuted PARAFAC models/Up_ling_func_2");
% permutedPARAFAC(microb_up_inter_cnt_scl,2,numReps,"./Permuted PARAFAC models/Up_inter_func_1");