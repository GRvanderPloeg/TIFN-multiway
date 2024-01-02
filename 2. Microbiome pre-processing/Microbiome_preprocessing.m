% PARAFAC functionality
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
% Split by niche
microb_tongue_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "tongue", :);
microb_lowling_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "lower jaw, lingual", :);
microb_lowinter_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "lower jaw, interproximal", :);
microb_upling_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "upper jaw, lingual", :);
microb_upinter_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "upper jaw, interproximal", :);
microb_saliva_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "saliva", :);

%%
% Select ASVs based on sparsity per RF group
sparsityThreshold = 50;

[tongue_sparsity_low, tongue_sparsity_mid, tongue_sparsity_high] = calculateRFsparsity(microb_tongue_raw, rf_low, rf_mid, rf_high);
[lowling_sparsity_low, lowling_sparsity_mid, lowling_sparsity_high] = calculateRFsparsity(microb_lowling_raw, rf_low, rf_mid, rf_high);
[lowinter_sparsity_low, lowinter_sparsity_mid, lowinter_sparsity_high] = calculateRFsparsity(microb_lowinter_raw, rf_low, rf_mid, rf_high);
[upling_sparsity_low, upling_sparsity_mid, upling_sparsity_high] = calculateRFsparsity(microb_upling_raw, rf_low, rf_mid, rf_high);
[upinter_sparsity_low, upinter_sparsity_mid, upinter_sparsity_high] = calculateRFsparsity(microb_upinter_raw, rf_low, rf_mid, rf_high);
[saliva_sparsity_low, saliva_sparsity_mid, saliva_sparsity_high] = calculateRFsparsity(microb_saliva_raw, rf_low, rf_mid, rf_high);

tongue_ASV_selection = (tongue_sparsity_low <= sparsityThreshold) | (tongue_sparsity_mid <= sparsityThreshold) | (tongue_sparsity_high <= sparsityThreshold);
lowling_ASV_selection = (lowling_sparsity_low <= sparsityThreshold) | (lowling_sparsity_mid <= sparsityThreshold) | (lowling_sparsity_high <= sparsityThreshold);
lowinter_ASV_selection = (lowinter_sparsity_low <= sparsityThreshold) | (lowinter_sparsity_mid <= sparsityThreshold) | (lowinter_sparsity_high <= sparsityThreshold);
upling_ASV_selection = (upling_sparsity_low <= sparsityThreshold) | (upling_sparsity_mid <= sparsityThreshold) | (upling_sparsity_high <= sparsityThreshold);
upinter_ASV_selection = (upinter_sparsity_low <= sparsityThreshold) | (upinter_sparsity_mid <= sparsityThreshold) | (upinter_sparsity_high <= sparsityThreshold);
saliva_ASV_selection = (saliva_sparsity_low <= sparsityThreshold) | (saliva_sparsity_mid <= sparsityThreshold) | (saliva_sparsity_high <= sparsityThreshold);

% Also remove nonbacterial ASVs
nonbacterial = ~((taxonomy(:,5) == "Chloroplast") | (taxonomy(:,6) == "Mitochondria"));

%%
% Separate metadata columns and numeric columns
microb_tongue_meta = microb_tongue_raw(:,1:5);
microb_tongue_numeric_strings = microb_tongue_raw(:,6:end);

microb_lowling_meta = microb_lowling_raw(:,1:5);
microb_lowling_numeric_strings = microb_lowling_raw(:,6:end);

microb_lowinter_meta = microb_lowinter_raw(:,1:5);
microb_lowinter_numeric_strings = microb_lowinter_raw(:,6:end);

microb_upling_meta = microb_upling_raw(:,1:5);
microb_upling_numeric_strings = microb_upling_raw(:,6:end);

microb_upinter_meta = microb_upinter_raw(:,1:5);
microb_upinter_numeric_strings = microb_upinter_raw(:,6:end);

microb_saliva_meta = microb_saliva_raw(:,1:5);
microb_saliva_numeric_strings = microb_saliva_raw(:,6:end);

%%
% Convert numeric data to actual numbers (so far they were strings)
microb_tongue_numeric = str2double(microb_tongue_numeric_strings);
microb_lowling_numeric = str2double(microb_lowling_numeric_strings);
microb_lowinter_numeric = str2double(microb_lowinter_numeric_strings);
microb_upling_numeric = str2double(microb_upling_numeric_strings);
microb_upinter_numeric = str2double(microb_upinter_numeric_strings);
microb_saliva_numeric = str2double(microb_saliva_numeric_strings);

%%
% CLR transformation
microb_tongue_clr = transformCLR(microb_tongue_numeric);
microb_lowling_clr = transformCLR(microb_lowling_numeric);
microb_lowinter_clr = transformCLR(microb_lowinter_numeric);
microb_upling_clr = transformCLR(microb_upling_numeric);
microb_upinter_clr = transformCLR(microb_upinter_numeric);
microb_saliva_clr = transformCLR(microb_saliva_numeric);

%%
% Reduce to the previously selected ASVs only
microb_tongue_reduced = microb_tongue_clr(:, (tongue_ASV_selection & nonbacterial));
microb_lowling_reduced = microb_lowling_clr(:, (lowling_ASV_selection & nonbacterial));
microb_lowinter_reduced = microb_lowinter_clr(:, (lowinter_ASV_selection & nonbacterial));
microb_upling_reduced = microb_upling_clr(:, (upling_ASV_selection & nonbacterial));
microb_upinter_reduced = microb_upinter_clr(:, (upinter_ASV_selection & nonbacterial));
microb_saliva_reduced = microb_saliva_clr(:, (saliva_ASV_selection & nonbacterial));

% Also reduce taxonomy to the same size
taxonomy_tongue_reduced = taxonomy((tongue_ASV_selection & nonbacterial), :);
taxonomy_lowling_reduced = taxonomy((lowling_ASV_selection & nonbacterial), :);
taxonomy_lowinter_reduced = taxonomy((lowinter_ASV_selection & nonbacterial), :);
taxonomy_upling_reduced = taxonomy((upling_ASV_selection & nonbacterial), :);
taxonomy_upinter_reduced = taxonomy((upinter_ASV_selection & nonbacterial), :);
taxonomy_saliva_reduced = taxonomy((saliva_ASV_selection & nonbacterial), :);

%%
% Reshape into a three-way matrix
%numTimepoints_microb = 7;
% microb_tongue = rawDataToCube_keepIndividuals(microb_tongue_reduced, microb_tongue_meta(:,2), microb_tongue_meta(:,3), numTimepoints_microb);
% microb_lowling = rawDataToCube_keepIndividuals(microb_lowling_reduced, microb_lowling_meta(:,2), microb_lowling_meta(:,3), numTimepoints_microb);
% microb_lowinter = rawDataToCube_keepIndividuals(microb_lowinter_reduced, microb_lowinter_meta(:,2), microb_lowinter_meta(:,3), numTimepoints_microb);
% microb_upling = rawDataToCube_keepIndividuals(microb_upling_reduced, microb_upling_meta(:,2), microb_upling_meta(:,3), numTimepoints_microb);
% microb_upinter = rawDataToCube_keepIndividuals(microb_upinter_reduced, microb_upinter_meta(:,2), microb_upinter_meta(:,3), numTimepoints_microb);
% microb_saliva = rawDataToCube_keepIndividuals(microb_saliva_reduced, microb_saliva_meta(:,2), microb_saliva_meta(:,3), numTimepoints_microb);

keepIndividuals = true;
microb_tongue = rawDataToCube(microb_tongue_reduced, microb_tongue_meta(:,2), str2double(microb_tongue_meta(:,3)), keepIndividuals);
microb_lowling = rawDataToCube(microb_lowling_reduced, microb_lowling_meta(:,2), str2double(microb_lowling_meta(:,3)), keepIndividuals);
microb_lowinter = rawDataToCube(microb_lowinter_reduced, microb_lowinter_meta(:,2), str2double(microb_lowinter_meta(:,3)), keepIndividuals);
microb_upling = rawDataToCube(microb_upling_reduced, microb_upling_meta(:,2), str2double(microb_upling_meta(:,3)), keepIndividuals);
microb_upinter = rawDataToCube(microb_upinter_reduced, microb_upinter_meta(:,2), str2double(microb_upinter_meta(:,3)), keepIndividuals);
microb_saliva = rawDataToCube(microb_saliva_reduced, microb_saliva_meta(:,2), str2double(microb_saliva_meta(:,3)), keepIndividuals);

%%
% Center the data
path_start = "./test_run/";

[microb_tongue_cnt, microb_tongue_means] = centerData(microb_tongue, 1);
[microb_lowling_cnt, microb_lowling_means] = centerData(microb_lowling, 1);
[microb_lowinter_cnt, microb_lowinter_means] = centerData(microb_lowinter, 1);
[microb_upling_cnt, microb_upling_means] = centerData(microb_upling, 1);
[microb_upinter_cnt, microb_upinter_means] = centerData(microb_upinter, 1);
[microb_saliva_cnt, microb_saliva_means] = centerData(microb_saliva, 1);

writematrix(microb_tongue_means, path_start+"Tongue_"+"means.csv");
writematrix(microb_lowling_means, path_start+"Lowling_"+"means.csv");
writematrix(microb_lowinter_means, path_start+"Lowinter_"+"means.csv");
writematrix(microb_upling_means, path_start+"Upling_"+"means.csv");
writematrix(microb_upinter_means, path_start+"Upinter_"+"means.csv");
writematrix(microb_saliva_means, path_start+"Saliva_"+"means.csv");

%%
% Scale the data
path_start = "./test_run/";

[microb_tongue_cnt_scl, microb_tongue_stds] = scaleData(microb_tongue_cnt, 2);
[microb_lowling_cnt_scl, microb_lowling_stds] = scaleData(microb_lowling_cnt, 2);
[microb_lowinter_cnt_scl, microb_lowinter_stds] = scaleData(microb_lowinter_cnt, 2);
[microb_upling_cnt_scl, microb_upling_stds] = scaleData(microb_upling_cnt, 2);
[microb_upinter_cnt_scl, microb_upinter_stds] = scaleData(microb_upinter_cnt, 2);
[microb_saliva_cnt_scl, microb_saliva_stds] = scaleData(microb_saliva_cnt, 2);

writematrix(microb_tongue_stds, path_start+"Tongue_"+"stds.csv");
writematrix(microb_lowling_stds, path_start+"Lowling_"+"stds.csv");
writematrix(microb_lowinter_stds, path_start+"Lowinter_"+"stds.csv");
writematrix(microb_upling_stds, path_start+"Upling_"+"stds.csv");
writematrix(microb_upinter_stds, path_start+"Upinter_"+"stds.csv");
writematrix(microb_saliva_stds, path_start+"Saliva_"+"stds.csv");

%%
% Create ID metadata
path_start = "./test_run/";

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
path_start = "./test_run/";

microb_tongue_ASV_meta = taxonomy_tongue_reduced(:, [2:end 1]);
microb_lowling_ASV_meta = taxonomy_lowling_reduced(:, [2:end 1]);
microb_lowinter_ASV_meta = taxonomy_lowinter_reduced(:, [2:end 1]);
microb_upling_ASV_meta = taxonomy_upling_reduced(:, [2:end 1]);
microb_upinter_ASV_meta = taxonomy_upinter_reduced(:, [2:end 1]);
microb_saliva_ASV_meta = taxonomy_saliva_reduced(:, [2:end 1]);

writematrix(microb_tongue_ASV_meta, path_start+"Tongue_"+"feature_meta.csv");
writematrix(microb_lowling_ASV_meta, path_start+"Lowling_"+"feature_meta.csv");
writematrix(microb_lowinter_ASV_meta, path_start+"Lowinter_"+"feature_meta.csv");
writematrix(microb_upling_ASV_meta, path_start+"Upling_"+"feature_meta.csv");
writematrix(microb_upinter_ASV_meta, path_start+"Upinter_"+"feature_meta.csv");
writematrix(microb_saliva_ASV_meta, path_start+"Saliva_"+"feature_meta.csv");

%%
% Save the cubes
path_start = "./test_run/";

writematrix(microb_tongue_cnt_scl, path_start+"Microb_tongue.csv");
writematrix(microb_lowling_cnt_scl, path_start+"Microb_lowling.csv");
writematrix(microb_lowinter_cnt_scl, path_start+"Microb_lowinter.csv");
writematrix(microb_upling_cnt_scl, path_start+"Microb_upling.csv");
writematrix(microb_upinter_cnt_scl, path_start+"Microb_upinter.csv");
writematrix(microb_saliva_cnt_scl, path_start+"Microb_saliva.csv");
