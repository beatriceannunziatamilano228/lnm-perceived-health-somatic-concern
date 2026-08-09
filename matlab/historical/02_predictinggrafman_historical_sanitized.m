% predictinggrafman.m
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

selectedColumns = [3:12 17:32 ,36];
lesion_size_result = SFArray(:, selectedColumns, end);


SF01_lesion_size = partialcorr(T0', SFArray(:,1),lesion_size_result); 
SF02_lesion_size = partialcorr(T0', SFArray(:,2),lesion_size_result);

NiiWrite2_NiiMask_local(SF01_lesion_size,'PROJECT_ROOT/3months/SF01map_312_1732 _lesion_size.nii',0,0);
NiiWrite2_NiiMask_local(SF02_lesion_size,'PROJECT_ROOT/3months/SF02map_312_1732 _lesion_size.nii',0,0);


unmask = zeros(902629,1);
unmask(find(NiiMask),:)=SF02_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
SF02_lesion_size_flipped = flip(unmask_3d, 1);
SF02_lesion_size_flipped = reshape(SF02_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(SF02_lesion_size_flipped,'PROJECT_ROOT/3months/SF02map_312_1732 _lesion_size_flipped.nii',0,1);

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



NBR01_lesion_size = partialcorr(T1', NBR(:,1), [NBR(:,2:end)], "Rows", "complete");
NBR02_lesion_size = partialcorr(T1', NBR(:,2), [NBR(:,[1 3:end])], "Rows", "complete");
NBR03_lesion_size = partialcorr(T1', NBR(:,3), [NBR(:,[1 2 4:end])], "Rows", "complete");
NBR04_lesion_size = partialcorr(T1', NBR(:,4), [NBR(:,[1 2 3 5:end])], "Rows", "complete");
NBR05_lesion_size = partialcorr(T1', NBR(:,5), [NBR(:,[1 2 3 4 6:end])], "Rows", "complete");
NBR06_lesion_size = partialcorr(T1', NBR(:,6), [NBR(:,[1 2 3 4 5 7:end])], "Rows", "complete");
NBR07_lesion_size = partialcorr(T1', NBR(:,7), [NBR(:,[1 2 3 4 5 6 8:end])], "Rows", "complete");
NBR08_lesion_size = partialcorr(T1', NBR(:,8), [NBR(:,[1 2 3 4 5 6 7 9:end])], "Rows", "complete");
NBR09_lesion_size = partialcorr(T1', NBR(:,9), [NBR(:,[1 2 3 4 5 6 7 8 10:end])], "Rows", "complete");
NBR10_lesion_size = partialcorr(T1', NBR(:,10), [NBR(:,[1 2 3 4 5 6 7 8 9 11:end])], "Rows", "complete");
NBR11_lesion_size = partialcorr(T1', NBR(:,11), [NBR(:,[1 2 3 4 5 6 7 8 9 10 12:end])], "Rows", "complete");
NBR12_lesion_size = partialcorr(T1', NBR(:,12), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 13:end])], "Rows", "complete");
NBR13_lesion_size = partialcorr(T1', NBR(:,13), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 14:end])], "Rows", "complete");
NBR14_lesion_size = partialcorr(T1', NBR(:,14), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 15:end])], "Rows", "complete");
NBR15_lesion_size = partialcorr(T1', NBR(:,15), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 16:end])], "Rows", "complete");
NBR16_lesion_size = partialcorr(T1', NBR(:,16), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17:end])], "Rows", "complete");
NBR17_lesion_size = partialcorr(T1', NBR(:,17), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18:end])], "Rows", "complete");
NBR18_lesion_size = partialcorr(T1', NBR(:,18), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19:end])], "Rows", "complete");
NBR19_lesion_size = partialcorr(T1', NBR(:,19), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20:end])], "Rows", "complete");
NBR20_lesion_size = partialcorr(T1', NBR(:,20), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21:end])], "Rows", "complete");
NBR21_lesion_size = partialcorr(T1', NBR(:,21), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22:end])], "Rows", "complete");
NBR22_lesion_size = partialcorr(T1', NBR(:,22), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23:end])], "Rows", "complete");
NBR23_lesion_size = partialcorr(T1', NBR(:,23), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24:end])], "Rows", "complete");
NBR24_lesion_size = partialcorr(T1', NBR(:,24), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25:end])], "Rows", "complete");
NBR25_lesion_size = partialcorr(T1', NBR(:,25), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26:end])], "Rows", "complete");
NBR26_lesion_size = partialcorr(T1', NBR(:,26), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27:end])], "Rows", "complete");
NBR27_lesion_size = partialcorr(T1', NBR(:,27), [NBR(:,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 28:end])], "Rows", "complete");


% NBR02
NiiWrite2_NiiMask_local(NBR02_lesion_size,'PROJECT_ROOT/Grafman/NBR02map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR02_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR02_lesion_size_flipped = flip(unmask_3d, 1);
NBR02_lesion_size_flipped = reshape(NBR02_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR02_lesion_size_flipped,'PROJECT_ROOT/3months/NBR02map_lesion_size_flipped.nii',0,1);
spat_r2 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR02_lesion_size_flipped));

% NBR03
NiiWrite2_NiiMask_local(NBR03_lesion_size,'PROJECT_ROOT/Grafman/NBR03map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR03_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR03_lesion_size_flipped = flip(unmask_3d, 1);
NBR03_lesion_size_flipped = reshape(NBR03_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR03_lesion_size_flipped,'PROJECT_ROOT/3months/NBR03map_lesion_size_flipped.nii',0,1);
spat_r3 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR03_lesion_size_flipped));

% NBR01
NiiWrite2_NiiMask_local(NBR01_lesion_size,'PROJECT_ROOT/Grafman/NBR01map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR01_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR01_lesion_size_flipped = flip(unmask_3d, 1);
NBR01_lesion_size_flipped = reshape(NBR01_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR01_lesion_size_flipped,'PROJECT_ROOT/3months/NBR01map_lesion_size_flipped.nii',0,1);
spat_r1 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR01_lesion_size_flipped));


% NBR04
NiiWrite2_NiiMask_local(NBR04_lesion_size,'PROJECT_ROOT/Grafman/NBR04map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR04_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR04_lesion_size_flipped = flip(unmask_3d, 1);
NBR04_lesion_size_flipped = reshape(NBR04_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR04_lesion_size_flipped,'PROJECT_ROOT/3months/NBR04map_lesion_size_flipped.nii',0,1);
spat_r4 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR04_lesion_size_flipped));

% NBR05
NiiWrite2_NiiMask_local(NBR05_lesion_size,'PROJECT_ROOT/Grafman/NBR05map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR05_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR05_lesion_size_flipped = flip(unmask_3d, 1);
NBR05_lesion_size_flipped = reshape(NBR05_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR05_lesion_size_flipped,'PROJECT_ROOT/3months/NBR05map_lesion_size_flipped.nii',0,1);
spat_r5 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR05_lesion_size_flipped));

% NBR06
NiiWrite2_NiiMask_local(NBR06_lesion_size,'PROJECT_ROOT/Grafman/NBR06map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR06_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR06_lesion_size_flipped = flip(unmask_3d, 1);
NBR06_lesion_size_flipped = reshape(NBR06_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR06_lesion_size_flipped,'PROJECT_ROOT/3months/NBR06map_lesion_size_flipped.nii',0,1);
spat_r6 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR06_lesion_size_flipped));

% NBR07
NiiWrite2_NiiMask_local(NBR07_lesion_size,'PROJECT_ROOT/Grafman/NBR07map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR07_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR07_lesion_size_flipped = flip(unmask_3d, 1);
NBR07_lesion_size_flipped = reshape(NBR07_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR07_lesion_size_flipped,'PROJECT_ROOT/3months/NBR07map_lesion_size_flipped.nii',0,1);
spat_r7 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR07_lesion_size_flipped));
% NBR08
NiiWrite2_NiiMask_local(NBR08_lesion_size,'PROJECT_ROOT/Grafman/NBR08map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR08_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR08_lesion_size_flipped = flip(unmask_3d, 1);
NBR08_lesion_size_flipped = reshape(NBR08_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR08_lesion_size_flipped,'PROJECT_ROOT/3months/NBR08map_lesion_size_flipped.nii',0,1);
spat_r8 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR08_lesion_size_flipped));
% NBR09
NiiWrite2_NiiMask_local(NBR09_lesion_size,'PROJECT_ROOT/Grafman/NBR09map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR09_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR09_lesion_size_flipped = flip(unmask_3d, 1);
NBR09_lesion_size_flipped = reshape(NBR09_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR09_lesion_size_flipped,'PROJECT_ROOT/3months/NBR09map_lesion_size_flipped.nii',0,1);
spat_r9 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR09_lesion_size_flipped));
% NBR10
NiiWrite2_NiiMask_local(NBR10_lesion_size,'PROJECT_ROOT/Grafman/NBR10map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR10_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR10_lesion_size_flipped = flip(unmask_3d, 1);
NBR10_lesion_size_flipped = reshape(NBR10_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR10_lesion_size_flipped,'PROJECT_ROOT/3months/NBR10map_lesion_size_flipped.nii',0,1);
spat_r10 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR10_lesion_size_flipped));
% NBR11
NiiWrite2_NiiMask_local(NBR11_lesion_size,'PROJECT_ROOT/Grafman/NBR11map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR11_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR11_lesion_size_flipped = flip(unmask_3d, 1);
NBR11_lesion_size_flipped = reshape(NBR11_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR11_lesion_size_flipped,'PROJECT_ROOT/3months/NBR11map_lesion_size_flipped.nii',0,1);
spat_r11 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR11_lesion_size_flipped));
% NBR12
NiiWrite2_NiiMask_local(NBR12_lesion_size,'PROJECT_ROOT/Grafman/NBR12map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR12_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR12_lesion_size_flipped = flip(unmask_3d, 1);
NBR12_lesion_size_flipped = reshape(NBR12_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR12_lesion_size_flipped,'PROJECT_ROOT/3months/NBR12map_lesion_size_flipped.nii',0,1);
spat_r12 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR12_lesion_size_flipped));
% NBR13
NiiWrite2_NiiMask_local(NBR13_lesion_size,'PROJECT_ROOT/Grafman/NBR13map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR13_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR13_lesion_size_flipped = flip(unmask_3d, 1);
NBR13_lesion_size_flipped = reshape(NBR13_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR13_lesion_size_flipped,'PROJECT_ROOT/3months/NBR13map_lesion_size_flipped.nii',0,1);
spat_r13 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR13_lesion_size_flipped));
% NBR14
NiiWrite2_NiiMask_local(NBR14_lesion_size,'PROJECT_ROOT/Grafman/NBR14map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR14_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR14_lesion_size_flipped = flip(unmask_3d, 1);
NBR14_lesion_size_flipped = reshape(NBR14_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR14_lesion_size_flipped,'PROJECT_ROOT/3months/NBR14map_lesion_size_flipped.nii',0,1);
spat_r14 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR14_lesion_size_flipped));
% NBR15
NiiWrite2_NiiMask_local(NBR15_lesion_size,'PROJECT_ROOT/Grafman/NBR15map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR15_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR15_lesion_size_flipped = flip(unmask_3d, 1);
NBR15_lesion_size_flipped = reshape(NBR15_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR15_lesion_size_flipped,'PROJECT_ROOT/3months/NBR15map_lesion_size_flipped.nii',0,1);
spat_r15 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR15_lesion_size_flipped));
% NBR16
NiiWrite2_NiiMask_local(NBR16_lesion_size,'PROJECT_ROOT/Grafman/NBR16map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR16_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR16_lesion_size_flipped = flip(unmask_3d, 1);
NBR16_lesion_size_flipped = reshape(NBR16_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR16_lesion_size_flipped,'PROJECT_ROOT/3months/NBR16map_lesion_size_flipped.nii',0,1);
spat_r16 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR16_lesion_size_flipped));
% NBR17
NiiWrite2_NiiMask_local(NBR17_lesion_size,'PROJECT_ROOT/Grafman/NBR17map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR17_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR17_lesion_size_flipped = flip(unmask_3d, 1);
NBR17_lesion_size_flipped = reshape(NBR17_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR17_lesion_size_flipped,'PROJECT_ROOT/3months/NBR17map_lesion_size_flipped.nii',0,1);
spat_r17 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR17_lesion_size_flipped));
% NBR18
NiiWrite2_NiiMask_local(NBR18_lesion_size,'PROJECT_ROOT/Grafman/NBR18map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR18_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR18_lesion_size_flipped = flip(unmask_3d, 1);
NBR18_lesion_size_flipped = reshape(NBR18_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR18_lesion_size_flipped,'PROJECT_ROOT/3months/NBR18map_lesion_size_flipped.nii',0,1);
spat_r18 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR18_lesion_size_flipped));
% NBR19
NiiWrite2_NiiMask_local(NBR19_lesion_size,'PROJECT_ROOT/Grafman/NBR19map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR19_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR19_lesion_size_flipped = flip(unmask_3d, 1);
NBR19_lesion_size_flipped = reshape(NBR19_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR19_lesion_size_flipped,'PROJECT_ROOT/3months/NBR19map_lesion_size_flipped.nii',0,1);
spat_r19 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR19_lesion_size_flipped));
% NBR20
NiiWrite2_NiiMask_local(NBR20_lesion_size,'PROJECT_ROOT/Grafman/NBR20map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR20_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR20_lesion_size_flipped = flip(unmask_3d, 1);
NBR20_lesion_size_flipped = reshape(NBR20_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR20_lesion_size_flipped,'PROJECT_ROOT/3months/NBR20map_lesion_size_flipped.nii',0,1);
spat_r20 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR20_lesion_size_flipped));
% NBR21
NiiWrite2_NiiMask_local(NBR21_lesion_size,'PROJECT_ROOT/Grafman/NBR21map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR21_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR21_lesion_size_flipped = flip(unmask_3d, 1);
NBR21_lesion_size_flipped = reshape(NBR21_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR21_lesion_size_flipped,'PROJECT_ROOT/3months/NBR21map_lesion_size_flipped.nii',0,1);
spat_r21 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR21_lesion_size_flipped));
% NBR22
NiiWrite2_NiiMask_local(NBR22_lesion_size,'PROJECT_ROOT/Grafman/NBR22map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR22_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR22_lesion_size_flipped = flip(unmask_3d, 1);
NBR22_lesion_size_flipped = reshape(NBR22_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR22_lesion_size_flipped,'PROJECT_ROOT/3months/NBR22map_lesion_size_flipped.nii',0,1);
spat_r22 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR22_lesion_size_flipped));
% NBR23
NiiWrite2_NiiMask_local(NBR23_lesion_size,'PROJECT_ROOT/Grafman/NBR23map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR23_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR23_lesion_size_flipped = flip(unmask_3d, 1);
NBR23_lesion_size_flipped = reshape(NBR23_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR23_lesion_size_flipped,'PROJECT_ROOT/3months/NBR23map_lesion_size_flipped.nii',0,1);
spat_r23 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR23_lesion_size_flipped));
% NBR24
NiiWrite2_NiiMask_local(NBR24_lesion_size,'PROJECT_ROOT/Grafman/NBR24map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR24_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR24_lesion_size_flipped = flip(unmask_3d, 1);
NBR24_lesion_size_flipped = reshape(NBR24_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR24_lesion_size_flipped,'PROJECT_ROOT/3months/NBR24map_lesion_size_flipped.nii',0,1);
spat_r24 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR24_lesion_size_flipped));
% NBR25
NiiWrite2_NiiMask_local(NBR25_lesion_size,'PROJECT_ROOT/Grafman/NBR25map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR25_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR25_lesion_size_flipped = flip(unmask_3d, 1);
NBR25_lesion_size_flipped = reshape(NBR25_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR25_lesion_size_flipped,'PROJECT_ROOT/3months/NBR25map_lesion_size_flipped.nii',0,1);
spat_r25 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR25_lesion_size_flipped));
% NBR26
NiiWrite2_NiiMask_local(NBR26_lesion_size,'PROJECT_ROOT/Grafman/NBR26map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR26_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR26_lesion_size_flipped = flip(unmask_3d, 1);
NBR26_lesion_size_flipped = reshape(NBR26_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR26_lesion_size_flipped,'PROJECT_ROOT/3months/NBR26map_lesion_size_flipped.nii',0,1);
spat_r26 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR26_lesion_size_flipped));
% NBR27
NiiWrite2_NiiMask_local(NBR27_lesion_size,'PROJECT_ROOT/Grafman/NBR27map_lesion_size.nii',0,0);
unmask = zeros(902629,1);
unmask(find(NiiMask),:)=NBR27_lesion_size;
unmask_3d = reshape(unmask,[91 109 91 1]);
NBR27_lesion_size_flipped = flip(unmask_3d, 1);
NBR27_lesion_size_flipped = reshape(NBR27_lesion_size_flipped, [], 1);
NiiWrite2_NiiMask_local(NBR27_lesion_size_flipped,'PROJECT_ROOT/3months/NBR27map_lesion_size_flipped.nii',0,1);
spat_r27 = corr(atanh(SF02_lesion_size_flipped), atanh(NBR27_lesion_size_flipped));

