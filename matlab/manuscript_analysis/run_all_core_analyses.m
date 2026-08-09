function run_all_core_analyses(dataDir, outputDir, nHistoricalPermutations, seed)
%RUN_ALL_CORE_ANALYSES Run the modular MATLAB analyses from restricted inputs.
% Peer-review Python analyses and SimNIBS simulations are run separately.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs"
    nHistoricalPermutations (1,1) double = 10000
    seed (1,1) double = 20260724
end

thisDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(thisDir));
referenceNifti = fullfile(repoRoot,"results","group_maps", ...
    "combined_lesion_associated_map_historical.nii.gz");
assert(isfile(referenceNifti),"Missing public reference NIfTI: %s",referenceNifti);
if ~isfolder(outputDir), mkdir(outputDir); end

run_primary_maps(dataDir,fullfile(outputDir,"primary_maps"),referenceNifti);
run_cross_dataset_correspondence(dataDir,fullfile(outputDir,"cross_dataset"), ...
    nHistoricalPermutations,seed);
run_cross_dataset_prediction(dataDir,fullfile(outputDir,"prediction"));
run_nrs_item_level_analysis(dataDir,fullfile(outputDir,"nrs_item_level"), ...
    referenceNifti);
run_historical_combined_map(dataDir,fullfile(outputDir,"combined_map"), ...
    referenceNifti);
run_supplementary_no_covariates(dataDir,fullfile(outputDir,"no_covariates"), ...
    referenceNifti);

controlFile = fullfile(dataDir,"Corbetta_additional_controls.mat");
if isfile(controlFile)
    run_dataset1_additional_controls(dataDir, ...
        fullfile(outputDir,"dataset1_controls"));
else
    warning("Optional additional-control file not found; skipping that analysis.");
end
end
