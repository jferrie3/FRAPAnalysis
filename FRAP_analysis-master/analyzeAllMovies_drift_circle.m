%% 
% put path to FRAP analysis code below
addpath C:\Users\User\Dropbox\Tjian_Darzacq\MATLAB_code\FRAP_analysis
% put path to BioFormats reader below
addpath c:/Users/User/Dropbox/MATLAB_code/bfmatlab/


%% circular bleach region

bleachFrame = 20;   % frame after which bleaching occurs
filterWidth = 5;    % width of Gaussian filter to use for smoothing
channelNum = 0;     % number of channel to use for defining FRAP ROI
radius = 5;
movieList = ls('../*.czi');
outdir = 'out';

mkdir(outdir);

for j = 1:size(movieList,1)
    fname = ['../' movieList(j,:)];
    currOutDir = [outdir '/' movieList(j,1:end-4)];
    mkdir(currOutDir)
    outfname = [currOutDir '/out'];
    pause(1);
    close all hidden;
    FRAPanalyze_drift_circle(fname,bleachFrame,filterWidth,channelNum,outfname,radius);
    
end