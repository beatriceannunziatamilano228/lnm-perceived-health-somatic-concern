function write_masked_nifti(vector, referenceNifti, outputFile)
%WRITE_MASKED_NIFTI Write a masked vector into the historical 2-mm MNI grid.
% The non-zero voxels of referenceNifti define the 285,903-voxel mask.

arguments
    vector double
    referenceNifti (1,1) string
    outputFile (1,1) string
end

info = niftiinfo(referenceNifti);
reference = niftiread(info);
mask = reference ~= 0;
assert(nnz(mask) == numel(vector), "Reference mask/vector length mismatch.");
volume = zeros(size(reference), "single");
volume(mask) = single(vector(:));

[outDir, base, ext] = fileparts(outputFile);
if ext == ".gz"
    [~, base2] = fileparts(base);
    base = base2;
end
if strlength(outDir) == 0, outDir = "."; end
if ~isfolder(outDir), mkdir(outDir); end

info.Datatype = "single";
info.BitsPerPixel = 32;
niftiwrite(volume, fullfile(outDir, base), info, "Compressed", true);
end
