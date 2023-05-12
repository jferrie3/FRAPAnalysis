% change this to wherever bioformats is on your computer
addpath F:/bfmatlab

%%
data = bfopen('F:/WT_Test_Set/20230110_Halo-WT_50nM-JF549_FRAP_cell001.czi');
% data = bfopen('20230116_Halo-WT_50nM-JF549_FRAP_cell018.czi');

%%
im = zeros(size(data{1}{1,1}));
for j = 1:size(data{1},1)
    im = im + double(data{1}{j,1});
end
    
%% run segmentation code
wNucleus = 5; % width of Gaussian filter for smoothing image to segment nucleus
wNucleoli = 1; % ditto for nucleoli
nuclearThresh = .5; % threshold for segmenting nucleus
nucleolarThresh = inf; % set this to infinity if you don't want to find any nucleoli
visualizationOn = 1; % set this to 0 if you want to run this in batch and don't want to see everything
[nm,nlm] = thresholdNucleusAndNucleoli(im,...
    wNucleus, wNucleoli,nuclearThresh,nucleolarThresh,visualizationOn)

%% use imfill if you want to fill in the holes in the mask
figure; imshow(imfill(nm,'holes'))


