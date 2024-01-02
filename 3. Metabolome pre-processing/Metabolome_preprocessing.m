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
% Load raw metabolomics data
metabolomics_raw = readmatrix("../0. Raw data input/20221005_wp2/metabolomics.tsv", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
metabolomicsHeader = metabolomics_raw(1,:);
metabolomics_raw_controlgroup = metabolomics_raw(ismember(metabolomics_raw(:,1), subjectsControl), :);

metabolomics_features = readmatrix("../0. Raw data input/20221005_wp2//metabolomics-features.tsv", Filetype="delimitedtext", Delimiter="\t", OutputType="string");

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
% Select metabolites based on their number of NAs
% Also consider the hand-made xenobiotics selection
na_threshold = 25;

metabolomics_na = metabolomics_raw_controlgroup(:,3:end) == "NA";
na_percentages =  sum(metabolomics_na) / size(metabolomics_na,1) * 100;
metabolomics_feature_selection_na = na_percentages <= na_threshold;

xenobiotics_data = readmatrix("../0. Raw data input/Xenobiotic_compounds_overview.csv", OutputType="string");
xenobiotics_selection = xenobiotics_data(xenobiotics_data(:,4) ~= "Remove", 2);
metabolomics_feature_selection_xeno = ismember(metabolomics_features(:,2), xenobiotics_selection);

metabolomics_feature_selection = metabolomics_feature_selection_na' | metabolomics_feature_selection_xeno;

%%
% Keep selected columns
metabolomics_reduced = metabolomics_raw_controlgroup(:, [true; true; metabolomics_feature_selection]);
metabolomics_features_reduced = metabolomics_features(metabolomics_feature_selection, :);

%%
% Separate metadata and numeric columns
% and convert numeric data to actual numbers (so far they were strings)
metabolomics_meta = metabolomics_reduced(:,1:2);
metabolomics_numeric_strings = metabolomics_reduced(:,3:end);
metabolomics_numeric = str2double(metabolomics_numeric_strings);

%%
% Impute the missing values with a random value between 0 and the detection
% limit
metabolomics_numeric_imp = metabolomics_numeric;
detectionLimits = min(metabolomics_numeric);
J = size(metabolomics_numeric, 2);

for j=1:J
    dataVector = metabolomics_numeric(:,j);
    missingData = isnan(dataVector);
    imputedData = rand(sum(missingData), 1) * detectionLimits(j);
    metabolomics_numeric_imp(missingData, j) = imputedData;
end

%%
% PQN
representativeSample = median(metabolomics_numeric_imp);
quotients = metabolomics_numeric_imp ./ representativeSample;
dilutionFactors = median(quotients,2);

halfDilutionFactors = dilutionFactors;
largerThanOne = dilutionFactors > 1;
smallerThanOne = dilutionFactors < 1;
halfDilutionFactors(largerThanOne) = halfDilutionFactors(largerThanOne) - ((halfDilutionFactors(largerThanOne) - 1) / 2);
halfDilutionFactors(smallerThanOne) = halfDilutionFactors(smallerThanOne) + (abs(halfDilutionFactors(smallerThanOne) - 1) / 2);

metabolomics_numeric_pqn = metabolomics_numeric_imp ./ dilutionFactors;

%%
% Log transform
metabolomics_numeric_log = log(metabolomics_numeric_pqn);

%%
% Reshape into a three-way matrix
numTimepoints_metab = 5;
metabolomics = rawDataToCube_keepIndividuals(metabolomics_numeric_log, metabolomics_meta(:,1), metabolomics_meta(:,2), numTimepoints_metab);

%%
% Center the data
[metabolomics_cnt, metabolomics_means] = centerData(metabolomics, 1);

%%
% Scale the data
[metabolomics_cnt_scl, metabolomics_stds] = scaleData(metabolomics_cnt, 2);

%%
% Create metadata
individual_meta = rf_data(ismember(rf_data(:,1), unique(metabolomics_meta(:,1))),1);
feature_meta = metabolomics_features_reduced(:, [3 4 2]);

%%
% Save the preprocessed data
path_start = "./";

writematrix(metabolomics_cnt_scl, path_start+"Metabolomics.csv");
writematrix(metabolomics_means, path_start+"Metabolomics_means.csv");
writematrix(metabolomics_stds, path_start+"Metabolomics_stds.csv");
writematrix(individual_meta, path_start+"Metabolomics_id_meta.csv");
writematrix(feature_meta, path_start+"Metabolomics_feature_meta.csv");