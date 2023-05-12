% % % % script for trying out FRAPanalyze

addpath c:/Users/User/Dropbox/MATLAB_code/bfmatlab/

%%
fname = '../FRAP_experiment_2019_07_25__10_50_41.czi';
bleachFrame = 10;
filterWidth = 5;
channelNum = 0;
outfname = 'test';

out = FRAPanalyze(fname,bleachFrame,filterWidth,channelNum,outfname);