function [Threshold, p_value_permuted, Thresholded_Brain_Behavior_Map] = Westfall_Young_Peak_Correction(Brain_Data, Behavioral_Data, N_Permutation)

%% Read me
% Current script identifies potential peaks in a brain-behavior correlation map
% It draws on the Westfall-Young correction method
% It takes in following input:
    % Brain data (n_voxel by n_participants)
    % Behavioral data (n_participants)
    % Number of permutations to run

% It yields following output:
    % Threshold, the top 95th percentile value from the permutated data
    % Permutation test p-value (where behavioral data is permutated)
    % Thresholded brain-behavior correlation map


mainDir = '<PROJECT_ROOT>/';addpath(fullfile(mainDir, 'Nifti_commands'));
addpath(fullfile(mainDir, 'NIfTI_20140122'));
addpath(fullfile(mainDir, 'CorbettaBehavior'));
addpath(fullfile(mainDir, '3months'));
addpath(fullfile(mainDir, 'Grafman'));
addpath(fullfile(mainDir, 'CorbettaBehavior/parfor_progress'));

load('<PROJECT_ROOT>/Grafman/Grafman.mat', 'Grafman2');
load('101new.mat');

Corbetta_Brain_Data = CorbettaFinal;
% Corbetta_Brain_Data = Corbetta_Brain_Data';
Brain_Data_1 = Corbetta_Brain_Data;
SF1_35 = readtable("Sf1_35final.xlsx");
SFArray = table2array(SF1_35);
selectedColumns = [3:12 17:32 ,36];
SelectedColumns_LS = SFArray(:, selectedColumns, end);
Behavioral_Data_0 = SelectedColumns_LS %covariates
Behavioral_Data_1 = SFArray(:,2); 

% Process Grafman Data
load('Grafman.mat');
T1 = Grafman2.maps;
Grafman_Brain_Data = T1;
NBR = Grafman2.NBR_items;
NBR = single(NBR);
 
Brain_Data_2 = Grafman_Brain_Data;
NBR = single([Grafman2.NBR_items]);
Grafman_Behavioral_Data = NBR(:, 2);
Behavioral_Data_2 = Grafman_Behavioral_Data;
Behavioral_Data_3 = [NBR(:,1) NBR(:,3:end)] %Cov

%% Calculate correlation and threshold
r1 = partialcorr(Brain_Data_1', Behavioral_Data_1, Behavioral_Data_0,  'rows', 'complete');
r2 = partialcorr(Brain_Data_2', Behavioral_Data_2,Behavioral_Data_3, 'rows', 'complete');

% Calculate sample size (i.e., weight)
weight1 = 101; % Sample size for Corbetta dataset
weight2 = 181; % Sample size for Grafman dataset

% Calculate the weighted average correlation
weighted_average_r = (r1 * weight1 + r2 * weight2) / (weight1 + weight2);


r_max = max(weighted_average_r);

%% Permutation Test
% Initialize an array to store permuted correlation coefficients
r_permuted_max = zeros(N_Permutation, 1);

% Perform the permutation test
for i = 1:N_Permutation
    % Set seed for reproducibility
    rng(i);

    % Randomly permute the Behavioral_Data for Corbetta dataset
    permuted_behavioral_data_1 = Behavioral_Data_1(randperm(length(Behavioral_Data_1)));
    % Calculate the partial correlation coefficient for the permuted Corbetta data
    r_permuted_1 = partialcorr(Brain_Data_1', permuted_behavioral_data_1, Behavioral_Data_0, 'rows', 'complete');

    % Randomly permute the Behavioral_Data for Grafman dataset
    permuted_behavioral_data_2 = Behavioral_Data_2(randperm(length(Behavioral_Data_2)));
    % Calculate the partial correlation coefficient for the permuted Grafman data
    r_permuted_2 = partialcorr(Brain_Data_2', permuted_behavioral_data_2, Behavioral_Data_3, 'rows', 'complete');

    % Calculate weighted average r for the permuted data
    weighted_average_r_permuted = (r_permuted_1 * weight1 + r_permuted_2 * weight2) / (weight1 + weight2);
    r_permuted_max(i) = max(weighted_average_r_permuted);

    fprintf('\b\b\b\b%3d%%', round(i/N_Permutation*100));
end

% Calculate the threshold for significance
Threshold = prctile(r_permuted_max, 95);
p_value_permuted = sum(r_max <= r_permuted_max)/N_Permutation;

%% Create a thresholded correlation matrix
Threshold_Index_Vector = (abs(weighted_average_r) >= abs(Threshold));
Thresholded_Map = zeros(size(weighted_average_r)); % Initialize a map with zeros

% Set values above or equal to the threshold to their actual values
Thresholded_Map(Threshold_Index_Vector) = weighted_average_r(Threshold_Index_Vector);

% Assign the brain behavior map
Thresholded_Brain_Behavior_Map = Thresholded_Map;
