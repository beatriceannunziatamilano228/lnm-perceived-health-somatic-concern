function data = load_lnm_inputs(dataDir)
%LOAD_LNM_INPUTS Load and validate the four restricted analysis matrices.

arguments
    dataDir (1,1) string
end

required = [ ...
    "Corbettaconnectivitymaps101.mat", ...
    "CorbettaSF02_controls.mat", ...
    "GrafmanmapsnoNaN.mat", ...
    "NBRwithoutNaN_NBR02column1.mat"];
for f = required
    assert(isfile(fullfile(dataDir, f)), "Missing restricted input: %s", f);
end

cMapsFile = load(fullfile(dataDir, required(1)), "CorbettaFinal");
cBehFile = load(fullfile(dataDir, required(2)), "question2_controls");
gMapsFile = load(fullfile(dataDir, required(3)), "Grafman2");
gBehFile = load(fullfile(dataDir, required(4)), "NBRwithoutNaN_NBR02column1");

data.dataset1.maps = double(cMapsFile.CorbettaFinal);
data.dataset1.behavior = double(cBehFile.question2_controls);
data.dataset2.maps = double(gMapsFile.Grafman2.mapsCopy);
data.dataset2.behavior = double(gBehFile.NBRwithoutNaN_NBR02column1);

assert(isequal(size(data.dataset1.maps), [285903, 101]));
assert(isequal(size(data.dataset1.behavior), [101, 28]));
assert(isequal(size(data.dataset2.maps), [285903, 181]));
assert(isequal(size(data.dataset2.behavior), [181, 28]));

% Dataset 1: column 1 SF-02; columns 2:27 SF-03–12 and SF-17–32;
% column 28 lesion size.
data.dataset1.outcome = data.dataset1.behavior(:,1);
data.dataset1.covariates = data.dataset1.behavior(:,2:end);

% Dataset 2 supplied order: NRS-02, NRS-01, NRS-03…NRS-27, lesion size.
data.dataset2.outcome = data.dataset2.behavior(:,1);
data.dataset2.covariates = data.dataset2.behavior(:,2:end);
data.dataset2.nrs = nan(181,27);
data.dataset2.nrs(:,2) = data.dataset2.behavior(:,1);
data.dataset2.nrs(:,1) = data.dataset2.behavior(:,2);
data.dataset2.nrs(:,3:27) = data.dataset2.behavior(:,3:27);
data.dataset2.lesionSize = data.dataset2.behavior(:,28);
end
