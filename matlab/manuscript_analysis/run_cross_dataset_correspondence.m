function summary = run_cross_dataset_correspondence(dataDir, outputDir, nPermutations, seed)
%RUN_CROSS_DATASET_CORRESPONDENCE Observed map similarity and historical permutation.
% The permutation reproduces the historical working logic: outcome and
% covariate rows are permuted together, independently within each dataset.
% The historical random seed was not preserved; the exact historical p=.04
% is therefore treated as a reported result, while the frozen reviewer test
% is supplied under python/reviewer_analyses/.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs/cross_dataset"
    nPermutations (1,1) double = 10000
    seed (1,1) double = 20260724
end
addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
data = load_lnm_inputs(dataDir);
if ~isfolder(outputDir), mkdir(outputDir); end

[r1, valid1] = partial_correlation_map(data.dataset1.maps, ...
    data.dataset1.outcome, data.dataset1.covariates);
[r2, valid2] = partial_correlation_map(data.dataset2.maps, ...
    data.dataset2.outcome, data.dataset2.covariates);
rawR = corr(r1, r2, "Rows", "complete");
[fisherR, nVoxels] = fisher_spatial_correlation(r1, r2);

rng(seed, "twister");
null = nan(nPermutations,1);
for i = 1:nPermutations
    p1 = randperm(size(data.dataset1.behavior,1));
    p2 = randperm(size(data.dataset2.behavior,1));
    a = partial_correlation_map(data.dataset1.maps, ...
        data.dataset1.outcome(p1), data.dataset1.covariates(p1,:));
    b = partial_correlation_map(data.dataset2.maps, ...
        data.dataset2.outcome(p2), data.dataset2.covariates(p2,:));
    null(i) = fisher_spatial_correlation(a,b);
end
historicalOneSidedP = (1 + nnz(null >= fisherR)) / (nPermutations + 1);

summary = table(nnz(valid1), nnz(valid2), nVoxels, rawR, fisherR, ...
    nPermutations, seed, historicalOneSidedP, ...
    'VariableNames', {"dataset1_n","dataset2_n","n_voxels", ...
    "spatial_r_raw","spatial_r_fisher_z","n_permutations","seed", ...
    "historical_style_one_sided_p"});
writetable(summary, fullfile(outputDir,"cross_dataset_summary.csv"));
save(fullfile(outputDir,"historical_style_null.mat"),"null","-v7.3");
end
