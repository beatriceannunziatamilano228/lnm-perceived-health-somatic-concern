function summary = run_cross_dataset_prediction(dataDir, outputDir)
%RUN_CROSS_DATASET_PREDICTION Dataset 1 map predicting Dataset 2 NRS-02.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs/prediction"
end
addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
data = load_lnm_inputs(dataDir);
if ~isfolder(outputDir), mkdir(outputDir); end

map1 = partial_correlation_map(data.dataset1.maps, ...
    data.dataset1.outcome, data.dataset1.covariates);
predictor = corr(map1, data.dataset2.maps, "Rows", "complete")';
valid = isfinite(predictor) & isfinite(data.dataset2.outcome) & abs(predictor)<1;
[rho,p] = corr(atanh(predictor(valid)), data.dataset2.outcome(valid), ...
    "Type", "Spearman", "Rows", "complete");
[pearsonR,pearsonP] = corr(atanh(predictor(valid)), data.dataset2.outcome(valid), ...
    "Type", "Pearson", "Rows", "complete");

summary = table(nnz(valid), rho, p, pearsonR, pearsonP, ...
    'VariableNames', {"analytic_n","spearman_rho","spearman_p_two_sided", ...
    "pearson_r","pearson_p_two_sided"});
writetable(summary, fullfile(outputDir,"figure3_prediction_summary.csv"));
save(fullfile(outputDir,"figure3_prediction_values.mat"),"predictor","valid");
end
