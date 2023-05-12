%Meagan FRAP plotting November 2020
%Combines functionality to call Thomas Graham's FRAPanalyze_drift_circle
%function for quantifying FRAP recovery, doing drift analysis, and normalization. 
%Then this code allows for grouping of data by condition, exponential fitting and plotting. 

clear all;
close all;
                   %%%%%%% Edit once %%%%%%%%
%Add paths to where Thomas' FRAP functions and bfmatlab are saved
addpath F:/bfmatlab
addpath F:/FRAP_analysis-master

outdir = 'F:/Full_FRAP_Dataset_v1/1451-A485';
output_name = 'results'; %specify the .mat file names of the FRAP output files from Thomas' code
%Folder names for all the outputs
mult_cond = {'1451-A485'}; %Plot recoveries and fits for multiple conditions to compare. 
    %'F:/Full_FRAP_Dataset_v1/1451', ...
    %'F:/Full_FRAP_Dataset_v1/1451-A485', ...
    %'F:/Full_FRAP_Dataset_v1/1451-WT', ...
    %'F:/Full_FRAP_Dataset_v1/1451-WT-A485', ...
    %'F:/Full_FRAP_Dataset_v1/Core', ...
    %'F:/Full_FRAP_Dataset_v1/Core-CTR', ...
    %'F:/Full_FRAP_Dataset_v1/CTR', ...
    %'F:/Full_FRAP_Dataset_v1/dCore', ...
    %'F:/Full_FRAP_Dataset_v1/H2B', ...
    %'F:/Full_FRAP_Dataset_v1/KI-A4', ...
    %'F:/Full_FRAP_Dataset_v1/KI-E7', ...
    %'F:/Full_FRAP_Dataset_v1/NLS', ...
    %'F:/Full_FRAP_Dataset_v1/NTR', ...
    %'F:/Full_FRAP_Dataset_v1/NTR-Core', ...
    %'F:/Full_FRAP_Dataset_v1/OLD_Core-CTR', ...
    %'F:/Full_FRAP_Dataset_v1/WT', ...
    %'F:/Full_FRAP_Dataset_v1/WT_3uM', ...
    %'F:/Full_FRAP_Dataset_v1/WT-A485'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  FOR THE USER TO EDIT PER ANALYSIS %%%%%%%%%%%%%%

ToDo = 1; %If ToDo = 0, analyze raw .czi files using Thomas' functions. 
          %If ToDo = 1, read in normalized FRAP data from output .mat
          %files, fit and make plots or single cells in a single condition. 
%All of the files/folders need to be with single quotes '' to work!!
myFolder = outdir; %This should contain your raw movies (.czi files)
myOutputs = outdir;
myExt = '/*_Halo-*'; %specify what subfolders to read in for plotting together, use Regex expression
myTxtCsvOutDir = append(outdir, '/');

compiled_Csv = 'F:/Full_FRAP_Dataset_v1/FRAP_Analysis_Output.csv';

bleachFrame = 21;
frameLength = 305.8; %exposure time in milliseconds 

%ANALYZING RAW MOVIES, Uses Thomas' FRAPanalyze_drift_circle function (ToDo = 0)
filterWidth = 10;    % width of Gaussian filter to use for smoothing
channelNum = 0;     % number of channel to use for defining FRAP ROI, starts at 0 for only 1 channel. 

%PLOTTING - plot single cells or fitted recoveries for a single or multiple conditions (ToDo = 1)
%single_cond = 'Torin1_nuc'; %Can be a single char or a char array. Uses Regex *mycondition* to select files of interest.
single_cell = 0; %=0 for all cells mean, =1 for single cell recoveries
single_cell_overlay = 1; %=0 to plot each single cell in its own plot, =1 to overlay all single cells in one plot. 
error_type = 0; %Error_type = 0 for standard deviation, =1 for 95% confidence intervals. NOTE: If you plot multiple conditions, it will ALWAYS plot standard deviation! 

%GROUPING BY CONDITION FOR PLOTTING
%Make sure single_cell = 0 (grouping only does conglomerate plotting)
%Error is always standard deviation. Only set either single_cond OR mult_cond! Use % to comment out the other.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You don't need to touch anything below here. :)


if ToDo ==0
    %%%%% Quantifying FRAP recovery from raw movies. 
    %
    for k = 1:length(outdir_list)
        outdir = outdir_list{k};
        myFolder = outdir; %This should contain your raw movies (.czi files)
        myOutputs = outdir;
        myExt = '/*_Halo-*'; %specify what subfolders to read in for plotting together, use Regex expression
        myTxtCsvOutDir = [outdir, '/WT_Single_Cell_FixBD_FixKa_W_'];
        movieList = dir(fullfile(myFolder, '*.czi'));
        if ~isfolder(outdir)
            mkdir(outdir);
        end
        
        for j = 1:size(movieList,1)
            fname = fullfile(myFolder, movieList(j).name);
            currOutDir = fullfile(outdir, movieList(j).name(:,1:end-4));
            if isfolder(currOutDir)
                continue
            end
            mkdir(currOutDir)
            outfname = fullfile(currOutDir, output_name);
            pause(1);
            close all hidden;
            %FRAPanalyze_drift_MNE(fname,bleachFrame,filterWidth,channelNum,outfname);
            FRAPanalyze_v1(fname,bleachFrame,filterWidth,channelNum,outfname,frameLength);

        end
        close all;
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting FRAP recoveries, fitting to exponential model, printing results
allMovies=[];
data = [];
data2=[];
movieFrames=[];

if ToDo ==1
    if exist('single_cond') ==1
        allMovies = dir(fullfile(myOutputs, ['*' single_cond '*'])); %loads files that match condition in filename
    else 
        mydir = dir(myOutputs);
        allMovies = mydir([mydir(:).isdir]);
        allMovies = allMovies(~ismember({allMovies(:).name},{'.','..'}));
        allMovies = allMovies(contains({allMovies.name}, mult_cond{1}));
        %loads all .mat files in subdirectories of myOutputs folder
    end 
    
    %Analyze SINGLE CELL FRAP recovery curves

    if single_cell ==1
        disp(['Plotting single cell FRAP recovery curves for ', num2str(length(allMovies)), ' cells.']);
        
        if single_cell_overlay ==1
           figure('position',[100 100 450 300]); hold on;
           title("Single cell recoveries", 'FontSize',12, 'FontName', 'Helvetica');
           ylabel('Normalized Fluorescence (AU)', 'FontSize',12, 'FontName', 'Helvetica');
           xlabel('Time post-bleach(seconds)', 'FontSize',12, 'FontName', 'Helvetica');
           legend;
           %axis([min(time) max(time) 0 1.2]);
           for i=1:length(allMovies)
                disp(['File:  ' allMovies(i).name]);
                data = load(fullfile(myOutputs, allMovies(i).name, [output_name '.mat']));
                movieFrames = length(data.out.fnorm);
                time = linspace(-bleachFrame*(frameLength/1000),((movieFrames-bleachFrame)*(frameLength/1000)), movieFrames);
                %Plot normalized FRAP curve
                norm_FRAP = data.out.fnorm;
                %figure('position',[100 100 450 300]); %[x y width height]
                plot(time, norm_FRAP, 'LineWidth', 3, 'DisplayName', ['Cell' allMovies(i).name(end-1:end)]);
                %legend(allMovies(i).name(end-1:end), 'FontSize',12, 'FontName', 'Helvetica');
           end
           hold off;
        else
            for i=1:length(allMovies)
                disp(['File:  ' allMovies(i).name]);
                data = load(fullfile(myOutputs, allMovies(i).name, [output_name '.mat']));
                movieFrames = length(data.out.fnorm);
                time = linspace(-bleachFrame*(frameLength/1000),((movieFrames-bleachFrame)*(frameLength/1000)), movieFrames);
                %Plot normalized FRAP curve
                norm_FRAP = data.out.fnorm;
                figure('position',[100 100 450 300]); %[x y width height]
                hold on;
                plot(time, norm_FRAP, 'k-', 'LineWidth', 3);
                title(["Cell  " allMovies(i).name(end-1:end)], 'FontSize',12, 'FontName', 'Helvetica');
                ylabel('Normalized Fluorescence (AU)', 'FontSize',12, 'FontName', 'Helvetica');
                xlabel('Time post-bleach(seconds)', 'FontSize',12, 'FontName', 'Helvetica');
                axis([min(time) max(time) 0 1.2]);
                hold off;
            end 
        end
     

    %CONGLOMERATE and analyze pooled cell FRAP recovery curves

    elseif single_cell==0
        disp(['Plotting mean FRAP recovery across ', num2str(length(allMovies)), ' cells.'])
        disp(['Files:  ' allMovies(:).name]);
        %CONGLOMERATE the data, fit to exponentials, and plot
        
        %%%% NEED TO FINISH

       if exist('mult_cond') ==1
            disp(["Sorting by: " convertCharsToStrings(mult_cond)]);
            figure('position',[100 100 500 500]); %[x y width height]
            title('Mean FRAP recovery comparison between conditions', 'FontSize',14, 'FontName', 'Helvetica');
            ylabel('Normalized Fluorescence (AU)', 'FontSize',14, 'FontName', 'Helvetica');
            xlabel('Time post-bleach(seconds)', 'FontSize',14, 'FontName', 'Helvetica');
            
            %colors
            cm=lines(length(mult_cond));
            hold all;
            
            for m=1:length(mult_cond)
                meanFRAP=[];
                stdFRAP=[];
                fNorm=[];
                movieFrames2 =[];
                data2=[];
                recovery_data=[];
                allMovies_cond=allMovies(contains({allMovies.name}, mult_cond{m}));

                for i=1:length(allMovies_cond)
                    data_load2 = load(fullfile(myOutputs, allMovies_cond(i).name, [output_name '.mat']));       
                    data2 = [data2  data_load2];
                    movieFrames2=[movieFrames2  length(data2(i).out.fnorm)];
                 end

                moviecrop=min(movieFrames2);
                data2=[];
                data_weights = [];

                for j=1:length(allMovies_cond)
                    data_load2 = load(fullfile(myOutputs, allMovies_cond(j).name, [output_name '.mat']));        
                    data2 = [data2  data_load2];        
                    fNorm = [fNorm  data2(j).out.fnorm(1:moviecrop)];
                    data_weights = [data_weights 1/var(data2(j).out.fnorm(1:bleachFrame))];
                end
                
                data_weights = data_weights./sum(data_weights);
                
                for k=1:moviecrop
                     %meanFRAP(k)= mean(fNorm(k,:));
                     meanFRAP(k)=sum(data_weights.*fNorm(k,:));
                     stdFRAP(k)=std(fNorm(k,:));                     
                end

                time2 = linspace(-bleachFrame*(frameLength/1000),((moviecrop-bleachFrame)*(frameLength/1000)), moviecrop);
                time2 = time2 - time2(bleachFrame);
                %time2 = linspace(0,(moviecrop*(frameLength/1000)), moviecrop); %Plots
                %time from beginning instead of bleach frame = 0
                time2_fit = time2(bleachFrame:end-1);
                meanFRAP_fit = meanFRAP(bleachFrame:end-1);
                FRAP_next_frame = meanFRAP(bleachFrame+1:end);
                Weights_fit = abs(meanFRAP_fit-FRAP_next_frame);
                %Weights_fit = 1 - (meanFRAP(bleachFrame:end-1) - min(meanFRAP(bleachFrame:end))) / (max(meanFRAP(bleachFrame:end)-min(meanFRAP(bleachFrame:end))));
                %Weights_fit(Weights_fit < 0.1) = 0.0;
                %disp(Weights_fit);

                %DO FITTING - single and double double-exponential

                f1 = fittype('1 - A*exp(-ka*x)');
                [OneExp_fit, OneExp_param] = fit(time2_fit', meanFRAP_fit', f1, 'Lower', [0.1 0.0001], 'Upper', [0.9 0.5], 'StartPoint', [0.35 0.025]); 
                OneExp_CI = confint(OneExp_fit);

                f2 = fittype('1 - A*exp(-ka*x) - (1 - B - A)*exp(-kb*x)');
                %Upper for 4th parameter was originally 0.05
                %[TwoExp_fit, TwoExp_param] = fit(time2_fit', meanFRAP_fit', f2, 'Lower', [0.005 0.005 0.005 0], 'Upper', [0.8 4 10 2], 'StartPoint', [0.5 1.1 0.05 0.0005]); 
                [TwoExp_fit, TwoExp_param] = fit(time2_fit', meanFRAP_fit', f2, 'Lower', [0.005 min(meanFRAP) 0.005 0], 'Upper', [0.8 min(meanFRAP) 10 2], 'StartPoint', [0.5 min(meanFRAP) 0.05 0.0005]); 
                TwoExp_CI = confint(TwoExp_fit);

                xFit1 = time2_fit;
                yFit1 = 1 -  OneExp_fit.A.* exp(-OneExp_fit.ka .* xFit1);
                yFit1Residual = yFit1 - meanFRAP_fit;
                xFit2 = time2_fit;
                yFit2 = 1 - TwoExp_fit.A .* exp(-TwoExp_fit.ka .* xFit2) - (1 - TwoExp_fit.B - TwoExp_fit.A) .* exp(-TwoExp_fit.kb .* xFit2);
                yFit2Residual = yFit2 - meanFRAP_fit;

                %Print out fit values
                Fit1_text(1) = {'1-Exp fit: 1-A*exp(-ka*t) '};
                Fit1_text(2) = {['A (mobile fraction) = ', num2str(OneExp_fit.A)]};
                Fit1_text(3) = {['ka = ', num2str(round(OneExp_fit.ka,2)), ' s^-1 or 1/k (residence time)= ', num2str(round(1/OneExp_fit.ka,2)), 's']};    
                %Fit1_text(4) = {['C (immobile fraction) = ', num2str(OneExp_fit.C)]};

                Fit2_text(1) = {'2-Exp fit: 1-A*exp(-ka*t)-(1-B-A) *exp(-kb*t)'};
                Fit2_text(2) = {['A = ', num2str(TwoExp_fit.A)]};
                Fit2_text(3) = {['ka = ', num2str(round(TwoExp_fit.ka,2)), ' s^-1 or 1/k (short residence time) = ', num2str(round(1/TwoExp_fit.ka,2)), 's']};
                Fit2_text(4) = {['B (Bleach Depth) = ', num2str(TwoExp_fit.B)]};
                Fit2_text(5) = {['kb = ', num2str(round(TwoExp_fit.kb,2)), ' s^-1 or 1/k (long residence time) = ', num2str(round(1/TwoExp_fit.kb,2)), 's']};
                %Fit2_text(6) = {['C (Immobile Fraction) = ', num2str(TwoExp_fit.C)]};

                %Plot mean FRAP recovery over time
                p1= errorbar(time2, meanFRAP, stdFRAP,'LineWidth', 1, 'Color', cm(m,:), 'DisplayName', 'Standard Error'); set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                p2 = plot(time2, meanFRAP, 'LineWidth', 3, 'DisplayName', ['Condition: ' mult_cond{m} ',  Num. cells = ' num2str(length(mult_cond{m}))], 'Color', cm(m,:));
                output_csv_name = append(myTxtCsvOutDir, 'FRAP_Output_Recovery_Curve.csv');
                csvwrite(output_csv_name, [time2' meanFRAP' stdFRAP']);

                %Plot the exponential fits
                p3 = plot(xFit1, yFit1, 'r--', 'LineWidth', 1, 'DisplayName', '1-Exp Fit'); set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                p4 = plot(xFit2, yFit2, 'g--', 'LineWidth', 1, 'DisplayName', '2-Exp Fit');  set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                output_csv_name = append(myTxtCsvOutDir, 'FRAP_Output_1exp_Fit.csv');
                csvwrite(output_csv_name, [xFit1' yFit1']);
                output_csv_name = append(myTxtCsvOutDir, 'FRAP_Output_2exp_Fit.csv');
                csvwrite(output_csv_name, [xFit2' yFit2']);
                %figure;
                %p5 = plot(xFit1, yFit1Residual, 'r-', 'LineWidth', 3, 'DisplayName', '1-Exp Fit Residual');
                %figure;
                %p6 = plot(xFit2, yFit2Residual, 'g-', 'LineWidth', 3, 'DisplayName', '2-Exp Fit Residual');
                output_csv_name = append(myTxtCsvOutDir, 'FRAP_Output_1exp_Residual.csv');
                csvwrite(output_csv_name, [xFit1' yFit1Residual']);
                output_csv_name = append(myTxtCsvOutDir, 'FRAP_Output_2exp_Residual.csv');
                csvwrite(output_csv_name, [xFit2' yFit2Residual']);

                % Grab the fit information
                fit_coeffs = coeffvalues(TwoExp_fit);
                fit_confints = confint(TwoExp_fit);

                write_params = [mult_cond{m} ',' num2str(fit_coeffs(1)) ',' num2str(fit_confints(1)) ',' num2str(fit_confints(2)) ',' num2str(fit_coeffs(2)) ',' num2str(fit_confints(3)) ',' num2str(fit_confints(4))  ',' num2str(fit_coeffs(3)) ',' num2str(fit_confints(5)) ',' num2str(fit_confints(6)) ',' num2str(fit_coeffs(4)) ',' num2str(fit_confints(7)) ',' num2str(fit_confints(8)) ',' num2str(TwoExp_param.rsquare) ',' num2str(TwoExp_param.rmse)];
                fileid = fopen(compiled_Csv, 'a');
                nbytes = fprintf(fileid, '%s\n', write_params);
                fclose(fileid);

                output_txt_name = append(myTxtCsvOutDir, 'FRAP_Fit_Info.txt');
                diary(output_txt_name);
                diary on;                
                disp(['----------------------------------------------------------'])
                disp(['Fitting Parameters and Values For ' [allMovies_cond(i).name]])
                disp(['----------------------------------------------------------'])
                disp([OneExp_fit])
                disp([OneExp_param])
                disp([TwoExp_fit])
                disp([TwoExp_param])
                diary off;

                %disp(['1-exp Long Res. time  ' [mult_cond{m}] '  ' num2str(round(1/OneExp_fit.ka,2))])            
                %disp(['2-exp Long Res. time  ' [mult_cond{m}] '  ' num2str(round(1/TwoExp_fit.kb,2))])            
                
                legend('Interpreter', 'none'); %This makes sure underscores aren't printed as subscripts          
                hold off;
            end
       else
            meanFRAP=[];
            fNorm=[];
            movieFrames2 =[];
            %Plot conglomerate recovery curves
            cm=hsv(length(allMovies)); %use "hsv" colomap to plot each cell as a different color
            figure('position',[100 100 500 500]); %[x y width height]
            hold on;

            for i=1:length(allMovies)
                data_load = load(fullfile(myOutputs, allMovies(i).name, [output_name '.mat']));       
                data = [data; data_load];        
                movieFrames2=[movieFrames2; length(data(i).out.fnorm)];
             end

            moviecrop=min(movieFrames2);

            for j=1:length(allMovies)
                data_load = load(fullfile(myOutputs, allMovies(i).name, [output_name '.mat']));        
                data = [data; data_load];        
                fNorm = [fNorm, data(j).out.fnorm(1:moviecrop)];
            end
            data=[];

            for k=1:moviecrop
                 meanFRAP(k)= mean(fNorm(k,:));
                 stdFRAP(k)=std(fNorm(k,:));
            end

            time2 = linspace(-bleachFrame*(frameLength/1000),((moviecrop-bleachFrame)*(frameLength/1000)), moviecrop);
            %time2 = linspace(0,(moviecrop*(frameLength/1000)), moviecrop); %Plots
            %time from beginning instead of bleach frame = 0

            time2_fit = time2(bleachFrame+1:end);
            meanFRAP_fit = meanFRAP(bleachFrame+1:end);

            %DO FITTING - single and double double-exponential

            f1 = fittype('1 - A*exp(-ka*x)');
            [OneExp_fit, OneExp_param] = fit(time2_fit', meanFRAP_fit', f1, 'Lower', [0.1 0.0001], 'Upper', [0.9 0.5], 'StartPoint', [0.35 0.025]); 
            OneExp_CI = confint(OneExp_fit);

            f2 = fittype('1 - A*exp(-ka*x) - B*exp(-kb*x)');
            [TwoExp_fit, TwoExp_param] = fit(time2_fit', meanFRAP_fit', f2, 'Lower', [0.05 0.05 0.05 0], 'Upper', [0.8 4 2 0.05], 'StartPoint', [0.5 1.1 0.5 0.0005]); 
            TwoExp_CI = confint(TwoExp_fit);

            xFit1 = min(time2):0.25:max(time2);
            yFit1 = 1 -  OneExp_fit.A.* exp(-OneExp_fit.ka .* xFit1);
            xFit2 = 0:0.25:max(time2);
            yFit2 = 1 - TwoExp_fit.A .* exp(-TwoExp_fit.ka .* xFit2) - TwoExp_fit.B .* exp(-TwoExp_fit.kb .* xFit2);

            %Print out fit values
            Fit1_text(1) = {'1-Exp fit: 1-A*exp(-ka*t) '};
            Fit1_text(2) = {['A (immobile fraction) = ', num2str(OneExp_fit.A)]};
            Fit1_text(3) = {['ka = ', num2str(round(OneExp_fit.ka,2)), ' s^-1 or 1/k (residence time)= ', num2str(round(1/OneExp_fit.ka,2)), 's']};    

            Fit2_text(1) = {'2-Exp fit: 1-A*exp(-ka*t)-B*exp(-kb*t)'};
            Fit2_text(2) = {['A = ', num2str(TwoExp_fit.A)]};
            Fit2_text(3) = {['ka = ', num2str(round(TwoExp_fit.ka,2)), ' s^-1 or 1/k (short residence time) = ', num2str(round(1/TwoExp_fit.ka,2)), 's']};
            Fit2_text(4) = {['B = ', num2str(TwoExp_fit.B)]};
            Fit2_text(5) = {['kb = ', num2str(round(TwoExp_fit.kb,2)), ' s^-1 or 1/k (long residence time) = ', num2str(round(1/TwoExp_fit.kb,2)), 's']};


            if error_type ==0 %Standard deviation error bars
                disp('Plot shows standard deviation of the mean.');
                hold on;
                %Plot mean FRAP recovery over time
                errorbar(time2, meanFRAP, stdFRAP,'o','LineWidth', 1);
                plot(time2, smooth(meanFRAP),'k-', 'LineWidth', 3);

                title(['Mean FRAP recovery of ', num2str(length(allMovies)), ' cells'], 'FontSize',14, 'FontName', 'Helvetica');
                ylabel('Normalized Fluorescence (AU)', 'FontSize',14, 'FontName', 'Helvetica');
                xlabel('Time post-bleach(seconds)', 'FontSize',14, 'FontName', 'Helvetica');
                axis([min(time2) max(time2) 0 1.2]);

                %Plot the fits
                plot(xFit1, yFit1, 'r-', 'LineWidth', 2);
                plot(xFit2, yFit2, 'g-', 'LineWidth', 2);

                %Overlay fit values as text on plot
                text(10,0.6,Fit1_text,'HorizontalAlignment','Left', 'FontSize',10, 'FontName', 'Helvetica');
                text(10,0.4,Fit2_text,'HorizontalAlignment','Left', 'FontSize',10, 'FontName', 'Helvetica');

                legend('Std. Dev', 'smoothed mean','1-Exp model fit', '2-Exp model fit','Location', 'SouthEast'); 
            end 


            if error_type ==1 %95% confidence interval error bars
                disp('Plot shows 95% confidence intervals of fits.')
                %Make exponential curves for the 95% confidence intervals
                yFit1_CI_l = 1 -  OneExp_CI(1,1).* exp(-OneExp_CI(1,2) .* xFit1);
                yFit1_CI_h = 1 -  OneExp_CI(2,1).* exp(-OneExp_CI(2,2) .* xFit1);


                yFit2_CI_l = 1 - TwoExp_CI(1,1) .* exp(-TwoExp_CI(1,3) .* xFit2) - TwoExp_CI(1,2) .* exp(-TwoExp_CI(1,4) .* xFit2);
                yFit2_CI_h = 1 - TwoExp_CI(2,1) .* exp(-TwoExp_CI(2,3) .* xFit2) - TwoExp_CI(2,2) .* exp(-TwoExp_CI(2,4) .* xFit2);

                %Plot mean FRAP recovery and overlay the fits
                hold on;
                %Plot mean FRAP recovery over time
                plot(time2, smooth(meanFRAP),'k-', 'LineWidth', 3);
                title(['Mean FRAP recovery of ', num2str(length(allMovies)), ' cells'], 'FontSize',14, 'FontName', 'Helvetica');
                ylabel('Normalized Fluorescence (AU)', 'FontSize',14, 'FontName', 'Helvetica');
                xlabel('Time post-bleach(seconds)', 'FontSize',14, 'FontName', 'Helvetica');
                axis([min(time2) max(time2) 0 1.2]);

                %Plot the fits
                plot(xFit1, yFit1, 'r-', 'LineWidth', 2);
                plot(xFit2, yFit2, 'g-', 'LineWidth', 2);

                %Plot the confidence interval lines
                plot(xFit1, yFit1_CI_l, 'r--', 'LineWidth', 1);
                plot(xFit1, yFit1_CI_h, 'r--', 'LineWidth', 1);

                plot(xFit2, yFit2_CI_l, 'g--', 'LineWidth', 1);
                plot(xFit2, yFit2_CI_h, 'g--', 'LineWidth', 1);

                %Shade the area between the 95% confidence intervals
                %1Exp fit
                xFit1_flip = [xFit1, fliplr(xFit1)];
                inBetween = [yFit1_CI_l, fliplr(yFit1_CI_h)];
                ci_1=fill(xFit1_flip, inBetween, 'r');
                set(ci_1,'facealpha',.1)
                %2Exp fit
                xFit2_flip = [xFit2, fliplr(xFit2)];
                inBetween = [yFit2_CI_l, fliplr(yFit2_CI_h)];
                ci_2=fill(xFit2_flip, inBetween, 'g');
                set(ci_2,'facealpha',.1)

                %Overlay fit values as text on plot
                text(10,0.3,Fit1_text,'HorizontalAlignment','Left', 'FontSize',9, 'FontName', 'Helvetica');
                text(10,0.6,Fit2_text,'HorizontalAlignment','Left', 'FontSize',9, 'FontName', 'Helvetica');

                legend('Smoothed mean', '1-Exp model fit', '2-Exp model fit','Location', 'SouthEast');
                hold off;

           end
       end
    end 
end

