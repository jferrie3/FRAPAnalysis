function [nuclearMask,nucleolarMask,allnucleiMask] = thresholdNucleusAndNucleoli(im,...
    wNucleus, wNucleoli,nuclearThresh,nucleolarThresh,visualizationOn)
% [nuclearMask,nucleolarMask] = thresholdNucleusAndNucleoli(im,...
%       wNucleus, wNucleoli, nuclearThresh,nucleolarThresh,visualizationOn)
% 
% outputs:
% nuclearMask - binary mask of nuclei
% nucleolarMask - binary mask of nucleoli
% 
% inputs:
% im - image (2d array)
% wNucleus - width of Gaussian filter for nuclear segmentation
% wNucleoli - width of Gaussian filter for nucleolar segmentation
% nuclearThresh - threshold to use for nuclear thresholding (number between
% 0 and 1, where 0 represents the minimum pixel intensity in the filtered
% image, and 1 represents the maximum intensity in the filtered image)
% nucleolarThresh - threshold to use for nucleolar thresholding 
% visualizationOn - set to 1 if you want to plot the output.
% 
% Thomas Graham, Tjian-Darzacq lab, 20200807

% apply a Gaussian filter (which can have a different width for detecting
% nuclei and detecting nucleoli
filteredForNucleus = gaussFilt2d(im,wNucleus);
filteredForNucleoli = gaussFilt2d(im,wNucleoli);

minPix = min(filteredForNucleus(:));
maxPix = mean(filteredForNucleus(:));
t1 = nuclearThresh * (maxPix - minPix) + minPix;

minPix = min(filteredForNucleoli(:));
maxPix = max(filteredForNucleoli(:));
t2 = nucleolarThresh* (maxPix - minPix) + minPix;

nuclearMask = filteredForNucleus > t1;
nucleolarMask = filteredForNucleoli > t2;

% Select only the largest nucleus in the frame
[labeledMask, numberofNuclei] = bwlabel(nuclearMask);
NucleiMeasurements = regionprops(labeledMask, 'area');
allAreas = [NucleiMeasurements.Area];
[sortedAreas, sortIndexed] = sort(allAreas, 'descend');
biggestNuclei = ismember(labeledMask, sortIndexed(1));
nuclearMask = biggestNuclei > 0;
allnucleiMask = labeledMask;

if visualizationOn
    figure;
    imagesc(im); hold on; colormap gray
    visboundaries(nuclearMask, 'Color','r')
    visboundaries(nucleolarMask,'Color','b')
    axis off;
    axis equal
end

end


function out = gaussFilt2d(im,sig)
% convolve an 2-dimensional image with a Gaussian filter
% usage: out = gaussFilt3d(im,sigXY,sigZ)
% sigXY and sigZ are standard deviation in the XY and Z directions

[y,x] = size(im);

[X,Y] = meshgrid(linspace(-(x-1)/2,(x-1)/2,x),linspace(-(y-1)/2,(y-1)/2,y));


g = exp(-(X.*X + Y.*Y)/(2*sig^2));
g = g/sum(g(:));

out = fftshift(ifft2(conj(fft2(g)).*fft2(im)));


end
