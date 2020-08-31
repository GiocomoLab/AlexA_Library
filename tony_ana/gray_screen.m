%%
sampRate = 4.1885;
%root = 'F:\LSD 2p first series\LSD 2p first series';
%sessionStructurePath = 'F:\LSD 2p first series\LSD 2p first series\Session_Structures';
%load('F:\LSD 2p first series\LSD 2p first series\Session_Structures\firstLSD_SessionNames.mat')
root = '/Users/attialex/Desktop/LSD 2p first series/';
sessionStructurePath = '/Users/attialex/Desktop/LSD 2p first series/Session_Structures';
load('/Users/attialex/Desktop/LSD 2p first series/Session_Structures/firstLSD_SessionNames.mat')
%% load data
GRAY_PRE={};
GRAY_POST = {};
BL_PRE={};
BL_POST={};

MAX_PRE={};
MAX_POST={};
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
        BL_PRE{end+1}=mean(act(:,1:S(1).greyStarts),2);
        MAX_PRE{end+1}=max(act(:,1:S(1).greyStarts),[],2);
        tmp =act(:,half_time+1+S(2).greyStarts(1):end);
        
        GRAY_POST{end+1} = tmp;
        start = half_time+1;
        stop = min(half_time+S(2).greyStarts(1),size(act,2));
        BL_POST{end+1}=mean(act(:,start:stop),2);
        MAX_POST{end+1}=max(act(:,start:stop),[],2);
        
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
        %rethrow(ME)
        warning(ME.message)
    end
end

%define analysis window
%% mean act

MPRE=[];
MPOST = [];
col = [];
for ii=1:length(GRAY_PRE)
    if size(GRAY_PRE{ii},2)>50 %throw out session with short baseline
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
end
figure
scatter(MPRE,MPOST,4,col)
ylabel('meanPost')
xlabel('meanPre')
axis equal
axis image
grid on
%% max act all
MPRE=[];
MPOST = [];
col = [];
cmap = jet(10);
colA=[];
for ii=1:length(GRAY_PRE)
    
    nc=size(GRAY_PRE{ii},1);
    MPRE = cat(1,MPRE,MAX_PRE{ii});
    MPOST = cat(1,MPOST,MAX_POST{ii});
    if LSD{ii}
        cc=[1 0 0];
    else
        cc=[0 0 1];
    end
    col = cat(1,col,repmat(cc,nc,1));
    colA=cat(1,colA,repmat(cmap(ii,:),nc,1));
    
end

figure
scatter(MPRE,MPOST,4,colA)
hold on
plot([0 100],[0 100],'k','LineWidth',2)
ylabel('maxPost')
xlabel('maxPre')
axis image
xlim([1 10])
ylim([1 10])

%% max act grey

MPRE=[];
MPOST = [];
col = [];
colA=[];
cmap = jet(10);
for ii=1:length(GRAY_PRE)
    if size(GRAY_POST{ii},2)>8 %throw out session with short baseline
    nc=size(GRAY_PRE{ii},1);
    MPRE = cat(1,MPRE,max(GRAY_PRE{ii},[],2));
    MPOST = cat(1,MPOST,max(GRAY_POST{ii},[],2));
    cc=LSD{ii};
    if LSD{ii}
        cc=[1 0 0];
    else
        cc=[0 0 1];
    end
    col = cat(1,col,repmat(cc,nc,1));
    colA=cat(1,colA,repmat(cmap(ii,:),nc,1));
    end
end
figure
scatter(MPRE,MPOST,4,colA)
hold on
plot([0 100],[0 100],'k','LineWidth',2)
ylabel('maxPost')
xlabel('maxPre')
axis image
xlim([1 10])
ylim([1 10])
grid on
%%

%% mean act normalized by mean activity pre grayscreen
tmp_PRE=[];
tmp_POST=[];
MPRE=[];
MPOST = [];
col = [];
for ii=1:length(GRAY_PRE)
    if size(GRAY_PRE{ii},2)>50
        tmp_PRE=cat(1,tmp_PRE,BL_PRE{ii});
        tmp_POST=cat(1,tmp_POST,BL_POST{ii});
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
end
figure

PRE = MPRE./tmp_PRE;
POST = MPOST./tmp_POST;
idx = randperm(numel(PRE));
scatter(PRE(idx),POST(idx),8,col(idx,:))
ylabel('meanPost')
xlabel('meanPre')
axis equal
axis image
grid on
box off

%% % time active
tmp_PRE=[];
tmp_POST=[];
MPRE=[];
MPOST = [];
col = [];
a_thresh = 1.1;
for ii=1:length(GRAY_PRE)
    if size(GRAY_PRE{ii},2)>50

        [nc,nf]=size(GRAY_PRE{ii});
        tmp_pre = sum(GRAY_PRE{ii}>a_thresh,2)/nf;
        tmp_post = sum(GRAY_POST{ii}>a_thresh,2)/nf;
        MPRE = cat(1,MPRE,tmp_pre);
        MPOST = cat(1,MPOST,tmp_post);
        [~,sidx]=sort(tmp_post,'descend');
        figure
        plot(GRAY_PRE{ii}(sidx(1),:))
        hold on
        plot(GRAY_POST{ii}(sidx(1),:))
        
        if LSD{ii}
            cc=[1 0 0];
        else
            cc=[0 0 1];
        end
        col = cat(1,col,repmat(cc,nc,1));
    end
end
figure


idx = randperm(numel(PRE));
scatter(MPRE,MPOST,8,col)
ylabel('meanPost')
xlabel('meanPre')
axis equal
axis image
grid on
box off


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
xlim([2 6])
ylim([2 6])
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
