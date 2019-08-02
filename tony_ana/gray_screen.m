%%
sampRate = 4.1885;
root = 'F:\LSD 2p first series\LSD 2p first series';
sessionStructurePath = 'F:\LSD 2p first series\LSD 2p first series\Session_Structures';
load('F:\LSD 2p first series\LSD 2p first series\Session_Structures\firstLSD_SessionNames.mat')
%% load data
GRAY_PRE={};
GRAY_POST = {};
GRAT_PRE = {};
SESSIONS = {};
LSD={};
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
    ACT = load(fullfile(root,exp,'actSmoothNorm'));  %actMinus07NeuPil; %actSmooth;
    tmp = fieldnames(ACT);
    act = ACT.(tmp{1});
    
    %figure
    %imagesc(act)
    half_time = size(act,2)/2;
    GRAY_PRE{end+1}= act(:,S(1).greyStarts(1):half_time);
    tmp =act(:,half_time+1+S(2).greyStarts(1):end);
    
    GRAY_POST{end+1} = tmp;
    SESSIONS{end+1}=expT;
    LSD{end+1}=SessionNames{expIdx,2};
    diffframes = abs((size(S(1).stimID,2)+ size(S(2).stimID,2))-size(act,2));
%     figure
%     plot(S(2).stimID')
%     hold on
%     plot(S(2).greyStarts(1),0,'ro')
    
    if diffframes>2
            warning(['frame number? diff: ' num2str(diffframes)])
    end
    catch ME
        warning(ME.message)
    end
end

%define analysis window
%% mean act
MPRE=[];
MPOST = [];
col = [];
for ii=1:length(GRAY_PRE)
    nc=size(GRAY_PRE{ii},1);
    MPRE = cat(1,MPRE,mean(GRAY_PRE{ii},2));
    MPOST = cat(1,MPOST,mean(GRAY_POST{ii},2));
    cc=LSD{ii};
    if LSD{ii}
        cc=[1 0 0];
    else
        cc=[0 0 1];
    end
    col = cat(1,col,repmat(cc,nc,1));
end
figure
scatter(MPRE,MPOST,4,col)
ylabel('meanPost')
xlabel('meanPre')
axis image
grid on
%% m
MPRE=[];
MPOST = [];
figure
for ii=1:length(GRAY_PRE)
    
    MPRE = mean(GRAY_PRE{ii},2);
    size(MPRE)
    MPOST = mean(GRAY_POST{ii},2);
    hold on
    if LSD{ii}
        col = 'r';
    else
        col = 'b';
    end
    errorbar(nanmean(MPRE),nanmean(MPOST),nanstd(MPRE),-nanstd(MPRE),nanstd(MPOST),-nanstd(MPOST),col)
end
ylabel('Post')
xlabel('Pre')
title('mean act')
axis square
axis image
grid on
%%
MPRE=[];
MPOST = [];
figure
for ii=1:length(GRAY_PRE)
    
    MPRE = max(GRAY_PRE{ii},[],2);
    size(MPRE)
    MPOST = max(GRAY_POST{ii},[],2);
    hold on
    if LSD{ii}
        col = 'r';
    else
        col = 'b';
    end
    if all(size(MPRE)==size(MPOST))
    errorbar(nanmean(MPRE),nanmean(MPOST),nanstd(MPRE),-nanstd(MPRE),nanstd(MPOST),-nanstd(MPOST),col)
    end
end
ylabel('Post')
xlabel('Pre')
title('max act')
grid on
axis image

%% pairwise correlations

MPRE = [];
MPOST = [];
figure
for ii=1:length(GRAY_PRE)
    nc=size(GRAY_PRE{ii},1);
    idx = triu(true(nc),1);
    MPRE = corrcoef(GRAY_PRE{ii}');
    MPRE=MPRE(idx);
    
    MPOST = corrcoef(GRAY_POST{ii}');
    MPOST = MPOST(idx);
    hold on
    if LSD{ii}
        col = 'r';
    else
        col = 'b';
    end
    if all(size(MPRE)==size(MPOST))
        errorbar(nanmean(MPRE),nanmean(MPOST),nanstd(MPRE),-nanstd(MPRE),nanstd(MPOST),-nanstd(MPOST),col)
    end
end
ylabel('Post')
xlabel('Pre')
title('pairwise correlation')
grid on
axis image
