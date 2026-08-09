function summary = run_supplementary_no_covariates(dataDir, outputDir, referenceNifti)
%RUN_SUPPLEMENTARY_NO_COVARIATES Compare primary maps with unadjusted maps.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs/no_covariates"
    referenceNifti (1,1) string = ""
end
addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
data = load_lnm_inputs(dataDir);
if ~isfolder(outputDir), mkdir(outputDir); end

r1Adjusted = partial_correlation_map(data.dataset1.maps, ...
    data.dataset1.outcome, data.dataset1.covariates);
r2Adjusted = partial_correlation_map(data.dataset2.maps, ...
    data.dataset2.outcome, data.dataset2.covariates);
r1Unadjusted = corr(data.dataset1.maps', data.dataset1.outcome, ...
    "Type","Pearson","Rows","complete");
r2Unadjusted = corr(data.dataset2.maps', data.dataset2.outcome, ...
    "Type","Pearson","Rows","complete");

r1 = corr(r1Adjusted,r1Unadjusted,"Rows","complete");
r2 = corr(r2Adjusted,r2Unadjusted,"Rows","complete");
summary = table(r1,r2,'VariableNames', ...
    {"dataset1_with_vs_without_covariates_spatial_r", ...
     "dataset2_with_vs_without_covariates_spatial_r"});
writetable(summary,fullfile(outputDir,"supplementary_S1_summary.csv"));
save(fullfile(outputDir,"no_covariate_maps.mat"), ...
    "r1Adjusted","r2Adjusted","r1Unadjusted","r2Unadjusted","-v7.3");

if strlength(referenceNifti)>0
    write_masked_nifti(r1Unadjusted, referenceNifti, ...
        fullfile(outputDir,"dataset1_no_covariates_r.nii.gz"));
    write_masked_nifti(r2Unadjusted, referenceNifti, ...
        fullfile(outputDir,"dataset2_no_covariates_r.nii.gz"));
end
end
