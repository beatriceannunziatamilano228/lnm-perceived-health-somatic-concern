function results = run_primary_maps(dataDir, outputDir, referenceNifti)
%RUN_PRIMARY_MAPS Reconstruct the two primary covariate-adjusted LNM maps.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs/primary_maps"
    referenceNifti (1,1) string = ""
end
addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
data = load_lnm_inputs(dataDir);
if ~isfolder(outputDir), mkdir(outputDir); end

[r1, valid1, rank1] = partial_correlation_map( ...
    data.dataset1.maps, data.dataset1.outcome, data.dataset1.covariates);
[r2, valid2, rank2] = partial_correlation_map( ...
    data.dataset2.maps, data.dataset2.outcome, data.dataset2.covariates);

results.dataset1_r = r1;
results.dataset2_r = r2;
results.dataset1_n = nnz(valid1);
results.dataset2_n = nnz(valid2);
results.dataset1_design_rank = rank1;
results.dataset2_design_rank = rank2;
save(fullfile(outputDir,"primary_maps.mat"),"-struct","results","-v7.3");

summary = table(results.dataset1_n, results.dataset2_n, rank1, rank2, ...
    'VariableNames', {"dataset1_complete_case_n","dataset2_complete_case_n", ...
    "dataset1_design_rank","dataset2_design_rank"});
writetable(summary, fullfile(outputDir,"primary_maps_summary.csv"));

if strlength(referenceNifti) > 0
    write_masked_nifti(r1, referenceNifti, ...
        fullfile(outputDir,"dataset1_primary_partial_r.nii.gz"));
    write_masked_nifti(r2, referenceNifti, ...
        fullfile(outputDir,"dataset2_primary_partial_r.nii.gz"));
end
end
