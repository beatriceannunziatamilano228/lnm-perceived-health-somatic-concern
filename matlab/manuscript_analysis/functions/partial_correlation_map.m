function [rMap, valid, designRank] = partial_correlation_map(maps, outcome, covariates)
%PARTIAL_CORRELATION_MAP Voxel-wise partial Pearson correlation.
% maps must be voxels x participants.

arguments
    maps double
    outcome double
    covariates double
end

outcome = outcome(:);
valid = isfinite(outcome) & all(isfinite(covariates),2);
assert(size(maps,2) == numel(outcome), "Participant dimension mismatch.");
assert(nnz(valid) > size(covariates,2) + 2, "Insufficient complete cases.");

designRank = rank([ones(nnz(valid),1), covariates(valid,:)]);
rMap = partialcorr(maps(:,valid)', outcome(valid), covariates(valid,:), ...
    "Type", "Pearson", "Rows", "complete");
end
