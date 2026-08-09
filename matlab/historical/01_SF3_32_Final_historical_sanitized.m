% SF3_32_Final.m
%
% RECOVERED HISTORICAL WORKING CODE — PROVENANCE ONLY
% Personal absolute paths were replaced with PROJECT_ROOT/ or USER_HOME/.
% Statistical operations were not intentionally modernized in this copy.
% The script depended on an interactive MATLAB workspace, local helper files,
% and restricted participant-level inputs; it is not the authoritative clean
% executable implementation. See ../clean/ for the modular version.
%

% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/Nifti_commands/
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/NIfTI_20140122/
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/CorbettaBehavior
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/3months
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/Grafman
load('PROJECT_ROOT/Nifti_commands/ROIs_old.mat','NiiMask');
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/CorbettaBehavior/parfor_progress/

% Make Corbetta Map

path = 'PROJECT_ROOT/Corbetta101';


files = dir(fullfile(path,'*.nii.gz'));
numFiles = length(files);
lesionSizes = zeros(numFiles,1);

for i = 1: numFiles
% Load the current lesion mask
mask = load_nii(fullfile(path,files(i).name));
% Count the number of positive voxels in the lesion mask
lesionSizes(i) = sum(mask.img(:) > 0);
end

load 101new.mat;

T0 = CorbettaFinal;
SF1_35 = readtable("Sf1_35final.xlsx");
SFArray = table2array(SF1_35);
selectedColumns = [3:32];
result = SFArray(:, selectedColumns, end);
selectedColumns = [3:32,36];
lesion_size_result = SFArray(:, selectedColumns, end);

% Mean centered covariates
mean_centered_result = result - mean(result,1);

% SF01 = partialcorr(T0', SFArray(:,1),result); %eliminating QS33-35
% SF02 = partialcorr(T0', SFArray(:,2),result);
% 
% NiiWrite2_NiiMask_local(SF01,'PROJECT_ROOT/3months/SF01map_3_32.nii',0,0);
% NiiWrite2_NiiMask_local(SF02,'PROJECT_ROOT/3months/SF02map_3_32.nii',0,0);

SF01_lesion_size = partialcorr(T0', SFArray(:,1),lesion_size_result); %eliminating QS33-36
SF02_lesion_size = partialcorr(T0', SFArray(:,2),lesion_size_result);
SF01_mean_centered = partialcorr(T0', SFArray(:,1) - mean(SFArray(:,1),1),mean_centered_result); %eliminating QS33-36
SF02_mean_centered = partialcorr(T0', SFArray(:,2) - mean(SFArray(:,2),1),mean_centered_result);

NiiWrite2_NiiMask_local(SF01_lesion_size,'PROJECT_ROOT/3months/SF01map_3_32_lesion_size.nii',0,0);
NiiWrite2_NiiMask_local(SF02_lesion_size,'PROJECT_ROOT/3months/SF02map_3_32_lesion_size.nii',0,0);
NiiWrite2_NiiMask_local(SF01_mean_centered,'PROJECT_ROOT/3months/SF01map_3_32_mean_centered.nii',0,0);
NiiWrite2_NiiMask_local(SF02_mean_centered,'PROJECT_ROOT/3months/SF02map_3_32_mean_centered.nii',0,0);


% Make Grafman Map

% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/Grafman/
path = 'PROJECT_ROOT/Grafman/lesions';

files = dir(fullfile(path,'*.nii'));
numFiles = length(files);
lesionSizes = zeros(numFiles,1);

for i = 1:numFiles
% Load the current lesion mask
mask = load_nii(fullfile(path,files(i).name));
% Count the number of positive voxels in the lesion mask
lesionSizes(i) = sum(mask.img(:) > 0);
end

load Grafman.mat

load('PROJECT_ROOT/Grafman/Grafman.mat', 'Grafman2');

T1 = Grafman2.maps;

lesionSizes = single(lesionSizes);

NBR = Grafman2.NBR_items;
NBR = single(NBR);
 
NBR = [NBR,lesionSizes];
NBR_denan = NBR;
NBR_denan(any(isnan(NBR_denan), 2), :) = [];
NBR_mean_centered = NBR - mean(NBR_denan,1);

NBR02 = partialcorr(T1', NBR(:,2), [NBR(:,1) NBR(:,3:end)], "Rows", "complete");
NBR02_mean_centered = partialcorr(T1', NBR_mean_centered(:,2), [NBR_mean_centered(:,1) NBR_mean_centered(:,3:end)], "Rows", "complete");

NiiWrite2_NiiMask_local(NBR02,'PROJECT_ROOT/Grafman/NBR02map.nii',0,0);
NiiWrite2_NiiMask_local(NBR02_mean_centered,'PROJECT_ROOT/Grafman/NBR02map_mean_centered.nii',0,0);

spat_r = corr(atanh(SF02), atanh(NBR02));
spat_r_mean_centered = corr(atanh(SF02_mean_centered), atanh(NBR02_mean_centered));





%72, 74, 40 (min max median)

% Spatial correlation of the maps of two datasets is 0.5697 vs 0.7223 

NBR02_LNM = load_nii('PROJECT_ROOT/Grafman/NBR02map.nii');
NBR02_LNM = NBR02_LNM.img;
NBR02_LNM = reshape(NBR02_LNM,[],1);
NBR02_LNM(NiiMask == 0) = [];
spat_r = corr(atanh(SF02), atanh(NBR02_LNM));

% Spatial Permutation 10000x 0.0451 vs 0.0061 ; 

s1 = size(SFArray, 1);
s2 = size(NBR, 1);

nperm = 10000;

shuffle = zeros(nperm, 1);
SF_prim = SFArray(:,2);
SF_cont = SFArray(:,[3:32 36]);
NBR_prim = NBR(:,2);
NBR_cont =[NBR(:,1) NBR(:,3:end)];

parfor_progress(nperm/5);
parfor i = 1:nperm
    rp1 = randperm(s1)';
    rp2 = randperm(s2)';
    a = partialcorr(T0', SF_prim(rp1), SF_cont(rp1), 'type', 'pearson', "rows", "complete");
    b = partialcorr(T1', NBR_prim(rp2), NBR_cont(rp2), 'type', 'pearson', "rows", "complete");
    shuffle(i) = corr(atanh(a), atanh(b), 'rows', 'complete');
    if mod(i,5) == 0
        parfor_progress;
    end
end
parfor_progress(0);

p = length(shuffle(shuffle>spat_r));
p = p/nperm


SF02pred = load_nii('SF02map_3_32_lesion_size.nii');
SF02pred = SF02pred.img;
SF02pred(NiiMask == 0) = [];
Predictor = corr(SF02pred', Grafman2.maps);
r = corr(atanh(Predictor)', Grafman2.NBR_items(:,2), 'type', 'Pearson', 'rows', 'complete');
%Fisher transf

[NBRpred,p] = corr(atanh(Predictor)',Grafman2.NBR_items, "rows","complete", "type", "Spearman");
scatter(atanh(Predictor)', Grafman2.NBR_items(:,2));

combinedmap = (atanh(SF02)*101+atanh(NBR02_LNM)*181)/282 %???
NiiWrite2_NiiMask_local(combinedmap,'PROJECT_ROOT/Grafman/combinedmap.nii',0,0);


%TFCE ANALYSIS


addpath 'PROJECT_ROOT/Tesi_Unipi'
addpath 'PROJECT_ROOT/CorbettaBehavior'
addpath 'PROJECT_ROOT/MatlabTFCE-master'
addpath 'PROJECT_ROOT/Tesi_Unipi'/CorbettaBehavior/parfor_progress/

load Nifti_commands/ROIs_old.mat
load 101new.mat

% load('SF36.mat', 'SF36CopyCopy') % no LS, no similar covariates 


CorbettaNifti = zeros([numel(NiiMask) size(CorbettaFinal, 2)]);
CorbettaNifti(NiiMask==1,:) = CorbettaFinal;

imgs = reshape(CorbettaNifti,[91,109,91, size(CorbettaFinal,2)]);
SF1_35 = readtable("Sf1_35final.xlsx");
SFArray = table2array(SF1_35);


selectedColumns = [3:32 36];
covariate = SFArray(:, selectedColumns, end);
% covariate = [ones(size(covariate,1)), covariate];
% covariate = SF36CopyCopy(:, 3:33);


[CorbettaSF02_TFCE1, CorbettaSF02_TFCE2] = matlab_tfce('regression',2,imgs,[],covariate,10000,2,0.5,26,.1,0);

min(min(min(CorbettaSF02_TFCE1{1,1}))) %0.0100 10000x
min(min(min(CorbettaSF02_TFCE2{1,1}))) 


%SPLIT-HALF CROSS-VALIDATION



%!!split half cross-validation for SF02!!
%                                                                                             
% PATH-SANITIZED HISTORICAL COMMAND: cd PROJECT_ROOT/3months
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/Nifti_commands/
% PATH-SANITIZED HISTORICAL COMMAND: addpath PROJECT_ROOT/Nifti_commands/NIfTI_20140122/
addpath 'PROJECT_ROOT/CorbettaBehavior'
load Nifti_commands/ROIs_old.mat
load('101new.mat')
path = 'PROJECT_ROOT/Corbetta101';
addpath 'PROJECT_ROOT/Tesi_Unipi'/CorbettaBehavior/parfor_progress/


T0 = table2array(CorbettaFinal); %connectivity maps 
SF1_35 = readtable("Sf1_35final.xlsx");
SFArray = table2array(SF1_35);
selectedColumns = [3:32];
result = SFArray(:, selectedColumns, end);

% SFArray(1, :) = []; %removed SF01 HERE!!
% 
% SFArray = [SFArray;lesionSizes];
result = transpose(result);

nperm = 10000;
pznumb = size(CorbettaFinal,2);
pzhalf = ceil(pznumb/2);

parfor_progress(nperm/5);

for i = 1:nperm
    rp = randperm(pznumb);
    randmaps = T0(:,rp);
    randSFArray = result(:,rp);  
    rmhalf1 = randmaps(:, 1:pzhalf);  %m stands for maps
    rahalf1 = randSFArray(:,1:pzhalf); %a stands for Array with answers
    rmhalf2 = randmaps(:, (pzhalf+1):end);
    rahalf2 = randSFArray(:, (pzhalf+1):end);
    SF02_half1 = partialcorr(rmhalf1', rahalf1(2,:)', rahalf1([3:end],:)'); 
    SF02_half2 = partialcorr(rmhalf2', rahalf2(2,:)', rahalf2([3:end],:)');
    [r(i), p(i)] = corr(atanh(SF02_half1), atanh(SF02_half2));

    if mod(i,5) == 0
       parfor_progress;
    end

    % pred1 = corr(SF02_half2, rmhalf1);
    % pred2 = corr(SF02_half1, rmhalf2);
    % pred = [pred1,pred2];
    
    % [r(i), p(i)] = partialcorr(pred', randSFArray(1,:)', randSFArray(2:end,:)');
end

parfor_progress(0);

histogram(r)

p_hist = 1-(sum(r>0)/nperm)

rmean = mean(r); %0.2843
pmean = mean(p); %0.0024

