function tableOut = run_dataset1_additional_controls(dataDir, outputDir)
%RUN_DATASET1_ADDITIONAL_CONTROLS Add disability, mood and sex covariates.
% For each control, the base and augmented maps are compared within the same
% complete-case subset, preventing missingness from driving map differences.

arguments
    dataDir (1,1) string
    outputDir (1,1) string = "outputs/dataset1_controls"
end
addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
data = load_lnm_inputs(dataDir);
controlFile = fullfile(dataDir,"Corbetta_additional_controls.mat");
assert(isfile(controlFile), "Missing restricted control file: %s", controlFile);
c = load(controlFile);
if ~isfolder(outputDir), mkdir(outputDir); end

spec = {
    "NIHSS", {"nih_total"};
    "FIM", {"fim_total"};
    "FAM", {"fam_tot"};
    "RNLI", {"rnltotal"};
    "GDS", {"gdss_score"};
    "Sex", {"gender_code"};
    "SIP domains", {"sip_body","sip_social","sip_mob","sip_com", ...
        "sip_emo","sip_house","sip_alert","sip_amb", ...
        "sip_psychosoc","sip_physical"}
    };

n = size(spec,1); analyticN = nan(n,1); nAdded = nan(n,1); spatialR = nan(n,1);
for i = 1:n
    fields = spec{i,2}; Z = nan(101,numel(fields));
    for j = 1:numel(fields)
        assert(isfield(c,fields{j}), "Missing control variable: %s", fields{j});
        Z(:,j) = c.(fields{j})(:);
    end
    available = all(isfinite(Z),2);
    baseCov = data.dataset1.covariates;
    [baseMap, validBase] = partial_correlation_map( ...
        data.dataset1.maps(:,available), data.dataset1.outcome(available), ...
        baseCov(available,:));
    [adjustedMap, validAdjusted] = partial_correlation_map( ...
        data.dataset1.maps(:,available), data.dataset1.outcome(available), ...
        [baseCov(available,:), Z(available,:)]);
    assert(nnz(validBase)==nnz(validAdjusted));
    analyticN(i)=nnz(validAdjusted); nAdded(i)=numel(fields);
    spatialR(i)=fisher_spatial_correlation(baseMap,adjustedMap);
end

tableOut = table(string(spec(:,1)),analyticN,nAdded,spatialR, ...
    'VariableNames', {"control","analytic_n","n_added_covariates", ...
    "spatial_r_base_vs_control_adjusted"});
writetable(tableOut,fullfile(outputDir,"dataset1_additional_controls.csv"));
end
