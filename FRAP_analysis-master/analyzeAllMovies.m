bleachFrame = 10;   % frame after which bleaching occurs
filterWidth = 5;    % width of Gaussian filter to use for smoothing
channelNum = 0;     % number of channel to use for defining FRAP ROI
outfname = 'test';
movieList = ls('../*.czi');
outdir = 'out';


mkdir(outdir);
for j = 1:size(movieList,1)
    fname = ['../' movieList(j,:)];
    currOutDir = [outdir '/' movieList(j,1:end-4)];
    mkdir(currOutDir)
    outfname = [currOutDir '/out.mat']; 
      
    pause(1);
    close all hidden;
    FRAPanalyze(fname,bleachFrame,filterWidth,channelNum,outfname);
    
end