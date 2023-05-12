allMovies = ls('out/Tai-mNG*');
allMovies = allMovies(3:end,:);

data = load(['out/' allMovies(1,:) '/out.mat']);

% % plot, fit to a single and double double-exponential, and save
allF = {};
allF2 = {};
mkdir out/fits

bleachFrame = 10;
ft = fittype('a + b*(1-exp(-x/k))');
ft2 = fittype('a + b0*(1-exp(-x/k0)) + b1*(1-exp(-x/k1))');
% coeffnames(ft)

for movieIndex = 1:size(allMovies,1)
    
    currMovie = allMovies(movieIndex,:);
    currMovieFile = ['../', currMovie, '.czi'];
    r = bfGetReader(currMovieFile);
    omeMeta = r.getMetadataStore();
    j = 0;
    c = [];
    t = [];
    
%     try

        while ~isempty(omeMeta.getPlaneDeltaT(0,j))
            t(end+1) = omeMeta.getPlaneDeltaT(0,j).value();
            j = j+1;
            try
                if isempty(omeMeta.getPlaneDeltaT(0,j))
                    break
                end
            catch
                break
            end
        end
        t = t - t(bleachFrame+1); % time relative to bleach

        data = load(['out/', currMovie, '/out.mat']);

        nframes = size(data.out.fnorm,1);
        
        close all;
        figure; plot(t(1:nframes),data.out.fnorm(:,1),'bo')
        hold on
        plot(t(1:nframes),data.out.fnorm(:,2),'r-')
%     catch
%         allF{movieIndex} = [];
%         allF2{movieIndex} = [];
%     end
%     
    try
        f = fit(transpose(t(bleachFrame+1:nframes)), data.out.fnorm(bleachFrame+1:end,1),...
            ft,'StartPoint',[0.5,0.9,20])
        allF{movieIndex} = f;
        plot(t(bleachFrame+1:nframes),f(t(bleachFrame+1:nframes)),'b-','LineWidth',2)
    catch
        allF{movieIndex} = [];
    end
    try
        f2 = fit(transpose(t(bleachFrame+1:nframes)), data.out.fnorm(bleachFrame+1:end,1),...
            ft2,'StartPoint',[0.6,0.2,30,0.2,10])
        allF2{movieIndex} = f2;
        plot(t(bleachFrame+1:nframes),f2(t(bleachFrame+1:nframes)),'k-','LineWidth',2)
    catch
        allF2{movieIndex} = [];
    end
 
        
    drawnow;
    saveas(gcf,['out/fits/' num2str(movieIndex) '.fig'])
    
end

save('out/fits/allFits.mat','allF','allF2')

%%
load('out/fits/allFits.mat')
notEmptyF = {allF{cellfun(@(x)~isempty(x),allF)}};
allA = cellfun(@(x)x.a,notEmptyF);
allB = cellfun(@(x)x.b,notEmptyF);
allK = cellfun(@(x)x.k,notEmptyF);

% fraction bleached = 1-a
% asymptotic recovery = a + b
% recovery over starting point = (a + b) - a = b
% fractional recovery = b/(1-a)

%%
fractionRecovery = allB./(1-allA)
fractionBleached = 1-allA

%%
% nonOutliers = 1:21;
% nonOutliers(14) = [];
% 
% noPonA = 1:9;
% withPonA = 11:19;

noPonA = [1,3:8];
withPonA = [11:15,18:19];


%%

subplot(1,3,1)
% bar(mean(fractionRecovery(nonOutliers)))
hold on;
% errorbar(mean(fractionRecovery(nonOutliers)),std(fractionRecovery(nonOutliers)),...
%     'LineWidth',2,'Color','k')
plot(ones(size(fractionBleached(noPonA))),...
    fractionBleached(noPonA),'.','MarkerSize',10,'Color','b')
plot(1+[-.1,.1],[1,1]*mean(fractionBleached(noPonA)),'LineWidth',3,'Color','k')

plot(2*ones(size(fractionBleached(withPonA))),...
    fractionBleached(withPonA),'.','MarkerSize',10,'Color','r')
plot(2+[-.1,.1],[1,1]*mean(fractionBleached(withPonA)),'LineWidth',3,'Color','k')

axis([0,3,0,1.3])
ylabel('Fraction bleached')


subplot(1,3,2)
% bar(mean(fractionRecovery(nonOutliers)))
hold on;
% errorbar(mean(fractionRecovery(nonOutliers)),std(fractionRecovery(nonOutliers)),...
%     'LineWidth',2,'Color','k')
plot(ones(size(fractionRecovery(noPonA))),...
    fractionRecovery(noPonA),'.','MarkerSize',10,'Color','b')
plot(1+[-.1,.1],[1,1]*mean(fractionRecovery(noPonA)),'LineWidth',3,'Color','k')

plot(2*ones(size(fractionRecovery(withPonA))),...
    fractionRecovery(withPonA),'.','MarkerSize',10,'Color','r')
plot(2+[-.1,.1],[1,1]*mean(fractionRecovery(withPonA)),'LineWidth',3,'Color','k')

axis([0,3,0,1.3])
ylabel('Fraction recovery')


subplot(1,3,3)
hold on;
% errorbar(mean(allK(nonOutliers)),std(allK(nonOutliers)),...
%     'LineWidth',2,'Color','k')
plot(ones(size(allK(noPonA))),...
    allK(noPonA),'.','MarkerSize',10,'Color','b')

plot(2*ones(size(allK(withPonA))),...
    allK(withPonA),'.','MarkerSize',10,'Color','r')

plot(1+[-.1,.1],[1,1]*mean(allK(noPonA)),'LineWidth',3,'Color','k')

plot(2+[-.1,.1],[1,1]*mean(allK(withPonA)),'LineWidth',3,'Color','k')


axis([0,3,0,100])
ylabel('Recovery time (s)')

set(gcf,'Position',[680   813   673   285])

legend('- PonA', '+ PonA')

saveas(gcf,'out/fits/fitComparison_badFitsRemoved.fig')

%%
ranksum(fractionBleached(noPonA),fractionBleached(withPonA))
ranksum(fractionRecovery(noPonA),fractionRecovery(withPonA))
ranksum(allK(noPonA),allK(withPonA))












