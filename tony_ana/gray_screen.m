%%
sampRate = 4.1885;
root = '/Users/attialex/Desktop/LSD 2p first series/';
sessionStructurePath = '/Users/attialex/Desktop/LSD 2p first series/Session_Structures';
load('/Users/attialex/Desktop/LSD 2p first series/Session_Structures/SessionNames.mat')
%% load data
GRAY_PRE={};
GRAY_POST = {};
SESSIONS = {};
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
    ACT = load(fullfile('/Users/attialex/Desktop/LSD 2p first series',exp,'actSmoothNorm'));  %actMinus07NeuPil; %actSmooth;
    tmp = fieldnames(ACT);
    act = ACT.(tmp{1});
    
    %figure
    %imagesc(act)
    half_time = size(act,2)/2;
    GRAY_PRE{end+1}= act(:,S(1).greyStarts(1):half_time);
    tmp =act(:,half_time+1+S(2).greyStarts(1):end);
    
    GRAY_POST{end+1} = tmp;
    SESSIONS{end+1}=expT;
    diffframes = abs((size(S(1).stimID,2)+ size(S(2).stimID,2))-size(act,2));
    figure
    plot(S(2).stimID')
    hold on
    plot(S(2).greyStarts(1),0,'ro')
    
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
for ii=1:length(GRAY_PRE)
    size(GRAY_PRE{ii})
    MPRE = cat(1,MPRE,mean(GRAY_PRE{ii},2));
    MPOST = cat(1,MPOST,mean(GRAY_POST{ii},2));
end
figure
plot(MPRE,MPOST,'.')
%% m
MPRE=[];
MPOST = [];
figure
for ii=1:length(GRAY_PRE)
    
    MPRE = mean(GRAY_PRE{ii},2);
    size(MPRE)
    MPOST = mean(GRAY_POST{ii},2);
    hold on
    
    errorbar(nanmean(MPRE),nanmean(MPOST),nanstd(MPRE),-nanstd(MPRE),nanstd(MPOST),-nanstd(MPOST))
end
%%
MPRE=[];
MPOST = [];
figure
for ii=1:length(GRAY_PRE)
    
    MPRE = max(GRAY_PRE{ii},[],2);
    size(MPRE)
    MPOST = max(GRAY_POST{ii},[],2);
    hold on
    if all(size(MPRE)==size(MPOST))
    errorbar(nanmean(MPRE),nanmean(MPOST),nanstd(MPRE),-nanstd(MPRE),nanstd(MPOST),-nanstd(MPOST))
    end
end
%%