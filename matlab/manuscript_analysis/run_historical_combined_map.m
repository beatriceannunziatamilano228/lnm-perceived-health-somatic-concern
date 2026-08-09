function results = run_historical_combined_map(dataDir, outputDir, referenceNifti)
%RUN_HISTORICAL_COMBINED_MAP Sample-size-weighted historical group map.
% The recovered master working script Fisher-transformed both component maps
% and used nominal weights 101 and 181. The weighted Fisher-z values were
% written directly as the historical combined map.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs/combined_map"
    referenceNifti (1,1) string = ""
end
addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
data = load_lnm_inputs(dataDir);
if ~isfolder(outputDir), mkdir(outputDir); end
r1 = partial_correlation_map(data.dataset1.maps, ...
    data.dataset1.outcome, data.dataset1.covariates);
r2 = partial_correlation_map(data.dataset2.maps, ...
    data.dataset2.outcome, data.dataset2.covariates);
combinedZ = (101*atanh(r1) + 181*atanh(r2)) / 282;
results.dataset1_r = r1; results.dataset2_r = r2;
results.combined_weighted_fisher_z = combinedZ;
save(fullfile(outputDir,"historical_combined_map.mat"),"-struct","results","-v7.3");
if strlength(referenceNifti)>0
    write_masked_nifti(combinedZ, referenceNifti, ...
        fullfile(outputDir,"combined_weighted_fisher_z.nii.gz"));
end
end
