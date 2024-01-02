% PARAFAC functionality
home = ".";
cd(home)
addpath("..\Matlab scripts\Scripts\"); % own scripts
addpath("..\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
% Load raw microbiome data
microbiome_raw = readmatrix("../0. Raw data input/20221005_wp2/count-table.tsv", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
microbiome_raw_controlgroup = microbiome_raw(microbiome_raw(:,4) == "control", :);

taxonomy = readmatrix("../0. Raw data input/20221005_wp2/taxonomic-classification.tsv", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
taxonomyHeader = taxonomy(1,:);
taxonomy = taxonomy(2:end,:);

subjectsControl = sortrows(unique(microbiome_raw(:,[2 4]), "rows"), 1);
subjectsControl = subjectsControl(subjectsControl(:,2) == "control",1);

%%
% Load preprocessed data from R pipeline
numTimepoints_microb = 7;
suffix = "_functional_analysis_raw";

% Tongue
microb_tongue = readmatrix("Tongue" + suffix + ".csv");
%microb_tongue = reshape(microb_tongue, size(microb_tongue,1), size(microb_tongue,2)/numTimepoints_microb, numTimepoints_microb);
microb_tongue_meta = readmatrix("Tongue" + suffix + "_subjects.csv", OutputType="string");
microb_tongue_feature_meta = readmatrix("Tongue" + suffix + "_features.csv", OutputType="string");

% Lower jaw, lingual
microb_lowling = readmatrix("Lowling" + suffix + ".csv");
%microb_lowling = reshape(microb_low_ling, size(microb_low_ling,1), size(microb_low_ling,2)/numTimepoints_microb, numTimepoints_microb);
microb_lowling_meta = readmatrix("Lowling" + suffix + "_subjects.csv", OutputType="string");
microb_lowling_feature_meta = readmatrix("Lowling" + suffix + "_features.csv", OutputType="string");

% Lower jaw, interproximal
microb_lowinter = readmatrix("Lowinter" + suffix + ".csv");
%microb_lowinter = reshape(microb_low_inter, size(microb_low_inter,1), size(microb_low_inter,2)/numTimepoints_microb, numTimepoints_microb);
microb_lowinter_meta = readmatrix("Lowinter" + suffix + "_subjects.csv", OutputType="string");
microb_lowinter_feature_meta = readmatrix("Lowinter" + suffix + "_features.csv", OutputType="string");

% Upper jaw, lingual
microb_upling = readmatrix("Upling" + suffix + ".csv");
%microb_upling = reshape(microb_up_ling, size(microb_up_ling,1), size(microb_up_ling,2)/numTimepoints_microb, numTimepoints_microb);
microb_upling_meta = readmatrix("Upling" + suffix + "_subjects.csv", OutputType="string");
microb_upling_feature_meta = readmatrix("Upling" + suffix + "_features.csv", OutputType="string");

% Upper jaw, interproximal
microb_upinter = readmatrix("Upinter" + suffix + ".csv");
%microb_upinter = reshape(microb_up_inter, size(microb_up_inter,1), size(microb_up_inter,2)/numTimepoints_microb, numTimepoints_microb);
microb_upinter_meta = readmatrix("Upinter" + suffix + "_subjects.csv", OutputType="string");
microb_upinter_feature_meta = readmatrix("Upinter" + suffix + "_features.csv", OutputType="string");

% Saliva
microb_saliva = readmatrix("Saliva" + suffix + ".csv");
%microb_saliva = reshape(microb_saliva, size(microb_saliva,1), size(microb_saliva,2)/numTimepoints_microb, numTimepoints_microb);
microb_saliva_meta = readmatrix("Saliva" + suffix + "_subjects.csv", OutputType="string");
microb_saliva_feature_meta = readmatrix("Saliva" + suffix + "_features.csv", OutputType="string");

%%
% Mapping
% Get the KO selection list
KOselection = readmatrix("./KOselection.csv", OutputType="string");
%KO2pathwayName = readmatrix("../Data/kegg_KO2pathwayName.csv", OutputType="string");
%KO2pathwayName = KO2pathwayName(2:end,2:end);

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
rf_data_control = rf_data(ismember(rf_data(:,1), subjectsControl), :);
rf_low = rf_data(rf_data(:,2) == "0", 1);
rf_mid = rf_data(rf_data(:,2) == "1", 1);
rf_high = rf_data(rf_data(:,2) == "2", 1);

%%
% Test
sparsityThreshold = 50;

[tongue_sparsity_low, tongue_sparsity_mid, tongue_sparsity_high] = calculateRFsparsity([microb_tongue_meta microb_tongue], rf_low, rf_mid, rf_high);
[lowling_sparsity_low, lowling_sparsity_mid, lowling_sparsity_high] = calculateRFsparsity([microb_lowling_meta microb_lowling], rf_low, rf_mid, rf_high);
[lowinter_sparsity_low, lowinter_sparsity_mid, lowinter_sparsity_high] = calculateRFsparsity([microb_lowinter_meta microb_lowinter], rf_low, rf_mid, rf_high);
[upling_sparsity_low, upling_sparsity_mid, upling_sparsity_high] = calculateRFsparsity([microb_upling_meta microb_upling], rf_low, rf_mid, rf_high);
[upinter_sparsity_low, upinter_sparsity_mid, upinter_sparsity_high] = calculateRFsparsity([microb_upinter_meta microb_upinter], rf_low, rf_mid, rf_high);
[saliva_sparsity_low, saliva_sparsity_mid, saliva_sparsity_high] = calculateRFsparsity([microb_saliva_meta microb_saliva], rf_low, rf_mid, rf_high);

tongue_sparsity_selection = (tongue_sparsity_low <= sparsityThreshold) | (tongue_sparsity_mid <= sparsityThreshold) | (tongue_sparsity_high <= sparsityThreshold);
lowling_sparsity_selection = (lowling_sparsity_low <= sparsityThreshold) | (lowling_sparsity_mid <= sparsityThreshold) | (lowling_sparsity_high <= sparsityThreshold);
lowinter_sparsity_selection = (lowinter_sparsity_low <= sparsityThreshold) | (lowinter_sparsity_mid <= sparsityThreshold) | (lowinter_sparsity_high <= sparsityThreshold);
upling_sparsity_selection = (upling_sparsity_low <= sparsityThreshold) | (upling_sparsity_mid <= sparsityThreshold) | (upling_sparsity_high <= sparsityThreshold);
upinter_sparsity_selection = (upinter_sparsity_low <= sparsityThreshold) | (upinter_sparsity_mid <= sparsityThreshold) | (upinter_sparsity_high <= sparsityThreshold);
saliva_sparsity_selection = (saliva_sparsity_low <= sparsityThreshold) | (saliva_sparsity_mid <= sparsityThreshold) | (saliva_sparsity_high <= sparsityThreshold);

%%
% Remove uninformative columns
sparsityThreshold = 100;
variationThreshold = 0.025;

[tongue_variation_selection, microb_tongue_featureSparsity, microb_tongue_featureVariation] = removeColumns(microb_tongue, sparsityThreshold, variationThreshold);
[lowling_variation_selection, microb_low_ling_featureSparsity, microb_low_ling_featureVariation] = removeColumns( microb_lowling, sparsityThreshold, variationThreshold);
[lowinter_variation_selection, microb_low_inter_featureSparsity, microb_low_inter_featureVariation] = removeColumns(microb_lowinter, sparsityThreshold, variationThreshold);
[upling_variation_selection, microb_up_ling_featureSparsity, microb_up_ling_featureVariation] = removeColumns(microb_upling, sparsityThreshold, variationThreshold);
[upinter_variation_selection, microb_up_inter_featureSparsity, microb_up_inter_featureVariation] = removeColumns(microb_upinter, sparsityThreshold, variationThreshold);
[saliva_variation_selection, microb_saliva_featureSparsity, microb_saliva_featureVariation] = removeColumns(microb_saliva, sparsityThreshold, variationThreshold);

%%
% Diagnostic plot for sparsity threshold setting
% subplot(2,3,1); histogram(microb_tongue_featureSparsity); title("Tongue");
% subplot(2,3,2); histogram(microb_low_ling_featureSparsity); title("Low ling");
% subplot(2,3,3); histogram(microb_low_inter_featureSparsity); title("Low inter");
% subplot(2,3,4); histogram(microb_up_ling_featureSparsity); title("Up ling");
% subplot(2,3,5); histogram(microb_up_inter_featureSparsity); title("Up inter");
% subplot(2,3,6); histogram(microb_saliva_featureSparsity); title("Saliva");

%%
% Diagnostic plot for variation threshold setting
% subplot(2,3,1); histogram(microb_tongue_featureVariation); title("Tongue"); set(gca, 'Yscale', 'log'); ylim([1e-1 10e4]);
% subplot(2,3,2); histogram(microb_low_ling_featureVariation); title("Low ling"); set(gca, 'Yscale', 'log'); ylim([1e-1 10e4]);
% subplot(2,3,3); histogram(microb_low_inter_featureVariation); title("Low inter"); set(gca, 'Yscale', 'log'); ylim([1e-1 10e4]);
% subplot(2,3,4); histogram(microb_up_ling_featureVariation); title("Up ling"); set(gca, 'Yscale', 'log'); ylim([1e-1 10e4]);
% subplot(2,3,5); histogram(microb_up_inter_featureVariation); title("Up inter"); set(gca, 'Yscale', 'log'); ylim([1e-1 10e4]);
% subplot(2,3,6); histogram(microb_saliva_featureVariation); title("Saliva"); set(gca, "Yscale", "log"); ylim([1e-1 10e4]);

%%
% Diagnostic plot for distribution of some features
%[I,~,K] = size(microb_tongue_reduced);

%subplot(2,3,1); histogram(reshape(microb_tongue_reduced(:,176,:), I*K, 1)); title(microb_tongue_ASV_meta_reduced(176,2));
%subplot(2,3,2); histogram(reshape(microb_tongue_reduced(:,3774,:), I*K, 1)); title(microb_tongue_ASV_meta_reduced(3774,2));
%subplot(2,3,3); histogram(reshape(microb_tongue_reduced(:,3394,:), I*K, 1)); title(microb_tongue_ASV_meta_reduced(3394,2));
%subplot(2,3,4); histogram(reshape(microb_tongue_reduced(:,4093,:), I*K, 1)); title(microb_tongue_ASV_meta_reduced(4093,2));
%subplot(2,3,5); histogram(reshape(microb_tongue_reduced(:,3788,:), I*K, 1)); title(microb_tongue_ASV_meta_reduced(3788,2));
%subplot(2,3,6); histogram(reshape(microb_tongue_reduced(:,941,:), I*K, 1)); title(microb_tongue_ASV_meta_reduced(941,2));

%%
% TEST: feature selection first
tongueSelection = tongue_sparsity_selection & tongue_variation_selection & ismember(microb_tongue_feature_meta(:,1), KOselection);
microb_tongue_reduced = microb_tongue(:, tongueSelection);
microb_tongue_feature_meta_reduced = microb_tongue_feature_meta(tongueSelection, :);

lowlingSelection = lowling_sparsity_selection & lowling_variation_selection & ismember(microb_lowling_feature_meta(:,1), KOselection);
microb_lowling_reduced = microb_lowling(:, lowlingSelection);
microb_lowling_feature_meta_reduced = microb_lowling_feature_meta(lowlingSelection, :);

lowinterSelection = lowinter_sparsity_selection & lowinter_variation_selection & ismember(microb_lowinter_feature_meta(:,1), KOselection);
microb_lowinter_reduced = microb_lowinter(:, lowinterSelection);
microb_lowinter_feature_meta_reduced = microb_lowinter_feature_meta(lowinterSelection, :);

uplingSelection = upling_sparsity_selection & upling_variation_selection & ismember(microb_upling_feature_meta(:,1), KOselection);
microb_upling_reduced = microb_upling(:, uplingSelection);
microb_upling_feature_meta_reduced = microb_upling_feature_meta(uplingSelection, :);

upinterSelection = upinter_sparsity_selection & upinter_variation_selection & ismember(microb_upinter_feature_meta(:,1), KOselection);
microb_upinter_reduced = microb_upinter(:, upinterSelection);
microb_upinter_feature_meta_reduced = microb_upinter_feature_meta(upinterSelection, :);

salivaSelection = saliva_sparsity_selection & saliva_variation_selection & ismember(microb_saliva_feature_meta(:,1), KOselection);
microb_saliva_reduced = microb_saliva(:, salivaSelection);
microb_saliva_feature_meta_reduced = microb_saliva_feature_meta(salivaSelection, :);

%%
% TEST: CLR after
microb_tongue_clr = transformCLR(microb_tongue_reduced);
microb_lowling_clr = transformCLR(microb_lowling_reduced);
microb_lowinter_clr = transformCLR(microb_lowinter_reduced);
microb_upling_clr = transformCLR(microb_upling_reduced);
microb_upinter_clr = transformCLR(microb_upinter_reduced);
microb_saliva_clr = transformCLR(microb_saliva_reduced);

%%
% CLR transformation
%microb_tongue_clr = transformCLR(microb_tongue);
%microb_lowling_clr = transformCLR(microb_lowling);
%microb_lowinter_clr = transformCLR(microb_lowinter);
%microb_upling_clr = transformCLR(microb_upling);
%microb_upinter_clr = transformCLR(microb_upinter);
%microb_saliva_clr = transformCLR(microb_saliva);

%%
% Feature selection
% tongueSelection = tongue_sparsity_selection & tongue_variation_selection & ismember(microb_tongue_feature_meta(:,1), KOselection);
% microb_tongue_reduced = microb_tongue_clr(:, tongueSelection);
% microb_tongue_feature_meta_reduced = microb_tongue_feature_meta(tongueSelection, :);
% 
% lowlingSelection = lowling_sparsity_selection & lowling_variation_selection & ismember(microb_lowling_feature_meta(:,1), KOselection);
% microb_lowling_reduced = microb_lowling_clr(:, lowlingSelection);
% microb_lowling_feature_meta_reduced = microb_lowling_feature_meta(lowlingSelection, :);
% 
% lowinterSelection = lowinter_sparsity_selection & lowinter_variation_selection & ismember(microb_lowinter_feature_meta(:,1), KOselection);
% microb_lowinter_reduced = microb_lowinter_clr(:, lowinterSelection);
% microb_lowinter_feature_meta_reduced = microb_lowinter_feature_meta(lowinterSelection, :);
% 
% uplingSelection = upling_sparsity_selection & upling_variation_selection & ismember(microb_upling_feature_meta(:,1), KOselection);
% microb_upling_reduced = microb_upling_clr(:, uplingSelection);
% microb_upling_feature_meta_reduced = microb_upling_feature_meta(uplingSelection, :);
% 
% upinterSelection = upinter_sparsity_selection & upinter_variation_selection & ismember(microb_upinter_feature_meta(:,1), KOselection);
% microb_upinter_reduced = microb_upinter_clr(:, upinterSelection);
% microb_upinter_feature_meta_reduced = microb_upinter_feature_meta(upinterSelection, :);
% 
% salivaSelection = saliva_sparsity_selection & saliva_variation_selection & ismember(microb_saliva_feature_meta(:,1), KOselection);
% microb_saliva_reduced = microb_saliva_clr(:, salivaSelection);
% microb_saliva_feature_meta_reduced = microb_saliva_feature_meta(salivaSelection, :);

%%
% Reshape into a three-way matrix
numTimepoints_microb = 7;

microb_tongue = rawDataToCube_keepIndividuals(microb_tongue_reduced, microb_tongue_meta(:,2), microb_tongue_meta(:,3), numTimepoints_microb);
microb_lowling = rawDataToCube_keepIndividuals(microb_lowling_reduced, microb_lowling_meta(:,2), microb_lowling_meta(:,3), numTimepoints_microb);
microb_lowinter = rawDataToCube_keepIndividuals(microb_lowinter_reduced, microb_lowinter_meta(:,2), microb_lowinter_meta(:,3), numTimepoints_microb);
microb_upling = rawDataToCube_keepIndividuals(microb_upling_reduced, microb_upling_meta(:,2), microb_upling_meta(:,3), numTimepoints_microb);
microb_upinter = rawDataToCube_keepIndividuals(microb_upinter_reduced, microb_upinter_meta(:,2), microb_upinter_meta(:,3), numTimepoints_microb);
microb_saliva = rawDataToCube_keepIndividuals(microb_saliva_reduced, microb_saliva_meta(:,2), microb_saliva_meta(:,3), numTimepoints_microb);

%%
% Center the data
path_start = "./20230618_run/";

[microb_tongue_cnt, microb_tongue_means] = centerData(microb_tongue, 1);
[microb_low_ling_cnt, microb_lowling_means] = centerData(microb_lowling, 1);
[microb_low_inter_cnt, microb_lowinter_means] = centerData(microb_lowinter, 1);
[microb_up_ling_cnt, microb_upling_means] = centerData(microb_upling, 1);
[microb_up_inter_cnt, microb_upinter_means] = centerData(microb_upinter, 1);
[microb_saliva_cnt, microb_saliva_means] = centerData(microb_saliva, 1);

writematrix(microb_tongue_means, path_start+"Tongue_"+"means.csv");
writematrix(microb_lowling_means, path_start+"Lowling_"+"means.csv");
writematrix(microb_lowinter_means, path_start+"Lowinter_"+"means.csv");
writematrix(microb_upling_means, path_start+"Upling_"+"means.csv");
writematrix(microb_upinter_means, path_start+"Upinter_"+"means.csv");
writematrix(microb_saliva_means, path_start+"Saliva_"+"means.csv");

%%
% Scale the data
path_start = "./20230618_run/";

[microb_tongue_cnt_scl, microb_tongue_stds] = scaleData(microb_tongue_cnt, 2);
[microb_lowling_cnt_scl, microb_lowling_stds] = scaleData(microb_low_ling_cnt, 2);
[microb_lowinter_cnt_scl, microb_lowinter_stds] = scaleData(microb_low_inter_cnt, 2);
[microb_upling_cnt_scl, microb_upling_stds] = scaleData(microb_up_ling_cnt, 2);
[microb_upinter_cnt_scl, microb_upinter_stds] = scaleData(microb_up_inter_cnt, 2);
[microb_saliva_cnt_scl, microb_saliva_stds] = scaleData(microb_saliva_cnt, 2);

writematrix(microb_tongue_stds, path_start+"Tongue_"+"stds.csv");
writematrix(microb_lowling_stds, path_start+"Lowling_"+"stds.csv");
writematrix(microb_lowinter_stds, path_start+"Lowinter_"+"stds.csv");
writematrix(microb_upling_stds, path_start+"Upling_"+"stds.csv");
writematrix(microb_upinter_stds, path_start+"Upinter_"+"stds.csv");
writematrix(microb_saliva_stds, path_start+"Saliva_"+"stds.csv");

%%
% Create ID metadata
path_start = "./20230618_run/";

microb_tongue_id_meta = rf_data_control(:,1);
microb_lowling_id_meta = rf_data_control(:,1);
microb_lowinter_id_meta = rf_data_control(:,1);
microb_upling_id_meta = rf_data_control(:,1);
microb_upinter_id_meta = rf_data_control(:,1);
microb_saliva_id_meta = rf_data_control(:,1);

writematrix(microb_tongue_id_meta, path_start+"Tongue_"+"id_meta.csv");
writematrix(microb_lowling_id_meta, path_start+"Lowling_"+"id_meta.csv");
writematrix(microb_lowinter_id_meta, path_start+"Lowinter_"+"id_meta.csv");
writematrix(microb_upling_id_meta, path_start+"Upling_"+"id_meta.csv");
writematrix(microb_upinter_id_meta, path_start+"Upinter_"+"id_meta.csv");
writematrix(microb_saliva_id_meta, path_start+"Saliva_"+"id_meta.csv");

%%
% Create feature metadata
path_start = "./20230618_run/";

writematrix(microb_tongue_feature_meta_reduced, path_start+"Tongue_"+"feature_meta.csv");
writematrix(microb_lowling_feature_meta_reduced, path_start+"Lowling_"+"feature_meta.csv");
writematrix(microb_lowinter_feature_meta_reduced, path_start+"Lowinter_"+"feature_meta.csv");
writematrix(microb_upling_feature_meta_reduced, path_start+"Upling_"+"feature_meta.csv");
writematrix(microb_upinter_feature_meta_reduced, path_start+"Upinter_"+"feature_meta.csv");
writematrix(microb_saliva_feature_meta_reduced, path_start+"Saliva_"+"feature_meta.csv");

%%
% Save the cubes
path_start = "./20230618_run/";

writematrix(microb_tongue_cnt_scl, path_start+"Microb_tongue_functional.csv");
writematrix(microb_lowling_cnt_scl, path_start+"Microb_lowling_functional.csv");
writematrix(microb_lowinter_cnt_scl, path_start+"Microb_lowinter_functional.csv");
writematrix(microb_upling_cnt_scl, path_start+"Microb_upling_functional.csv");
writematrix(microb_upinter_cnt_scl, path_start+"Microb_upinter_functional.csv");
writematrix(microb_saliva_cnt_scl, path_start+"Microb_saliva_functional.csv");
