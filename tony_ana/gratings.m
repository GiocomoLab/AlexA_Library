
sampRate = 4.1885;
root = 'F:\LSD 2p first series\LSD 2p first series';
sessionStructurePath = 'F:\LSD 2p first series\LSD 2p first series\Session_Structures';
load('F:\LSD 2p first series\LSD 2p first series\Session_Structures\firstLSD_SessionNames.mat')
%% load data
GRAT_PRE={};
GRAT_POST = {};

SESSIONS = {};
LSD={};
win = round([-3*sampRate 6*sampRate]);
for expIdx = 1:length(SessionNames) % for each session structure
    
    % get date, mouse, exp name
    expT =  SessionNames{expIdx};
    P = strfind(expT,'.');
    date = expT(1:8);
    mouse = expT(10:P-1);
    exp = [date,mouse];
    
    %Load Session Structure S and calcium activity act incl. planeIdx
    try
    load(fullfile(sessionStructurePath,[date,'_',mouse]));
    ACT = load(fullfile(root,exp,'actMinus07NeuPil'));  %actMinus07NeuPil; %actSmooth;
    tmp = fieldnames(ACT);
    act = ACT.(tmp{1});
    
    %figure
    %imagesc(act)
    half_time = size(act,2)/2;
    for iGrat = 1:4
    GRAT_PRE{expIdx,iGrat}=extract_snps(act,S(1).gratingStarts{1,iGrat},'win',win);
    GRAT_POST{expIdx,iGrat}=extract_snps(act,S(2).gratingStarts{1,iGrat}+half_time+1,'win',win);
    end
    catch ME
        warning(ME.message)
    end
end
%%
masterTime = linspace(win(1), win(2),size(GRAT_PRE{end},2))/sampRate;
params.winIDX=-100:100;
params.masterTime=masterTime;
params.xLim=[-3 6];

figure
for iGrat = 1:4
avgGratPre = [];
avgGratPost = [];
for ii=1:length(GRAT_PRE)
    avgGratPre = cat(1,avgGratPre,mean(GRAT_PRE{ii,iGrat},3));
    avgGratPost = cat(1,avgGratPost,mean(GRAT_POST{ii,iGrat},3));
end
subplot(1,4,iGrat)
plotAVGSEM(avgGratPre',gca,'parameters',params,'ms',1,'baseline',find(masterTime<0))
plotAVGSEM(avgGratPost',gca,'parameters',params,'ms',1,'baseline',find(masterTime<0),'col',[1 0 0])
end
%%
masterTime = linspace(win(1), win(2),size(GRAT_PRE{end},2))/sampRate;
params.winIDX=-100:100;
params.masterTime=masterTime;
params.xLim=[-3 6];
for ii=1:length(GRAT_PRE)
figure
for iGrat = 1:4

    avgGratPre = mean(GRAT_PRE{ii,iGrat},3);
    avgGratPost = mean(GRAT_POST{ii,iGrat},3);

subplot(1,4,iGrat)
plotAVGSEM(avgGratPre',gca,'parameters',params,'ms',1,'baseline',find(masterTime<0))
plotAVGSEM(avgGratPost',gca,'parameters',params,'ms',1,'baseline',find(masterTime<0),'col',[1 0 0])
end
end


