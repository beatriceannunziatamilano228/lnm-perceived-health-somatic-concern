function [r, nVoxels] = fisher_spatial_correlation(mapA, mapB)
%FISHER_SPATIAL_CORRELATION Pearson spatial correlation after Fisher z.

mapA = mapA(:); mapB = mapB(:);
valid = isfinite(mapA) & isfinite(mapB) & abs(mapA) < 1 & abs(mapB) < 1;
nVoxels = nnz(valid);
r = corr(atanh(mapA(valid)), atanh(mapB(valid)), ...
    "Type", "Pearson", "Rows", "complete");
end
