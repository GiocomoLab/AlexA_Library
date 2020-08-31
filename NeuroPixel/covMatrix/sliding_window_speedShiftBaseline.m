matfiles = dir('/Users/attialex/slidingWindow_baseline_100cm/*.mat');

regions = {'MEC','VISp','RS'};
PEAKS = {};
SHIFTS = {};
SPEED = [];
REG={};
CLU=[];
FID = [];
SESSION = {};
for iF=1:numel(matfiles)
    
    [~,sn]=fileparts(matfiles(iF).name);
    parts = strsplit(sn,'_');
    
    
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    nC=size(data_out.corrMat,1);
    stab = nanmean(nanmean(data_out.corrMat,2),3);
    reg = data_out.region;
    idx = startsWith(reg,'VISpm');
    reg(idx)={'VISm'};
    idx = startsWith(reg,'RSPag');
    reg(idx)={'aglRS'};
    idx = startsWith(reg,'RS');
    reg(idx)={'RSC'};
    
    valid_idx = stab'>=0.5 & startsWith(reg,regions);
    if nnz(valid_idx)<5
        continue
    end
    parts = strsplit(sn,'_');
    reg=reg(valid_idx);
    %PEAKS = cat(1,PEAKS,data_out.PEAKS(valid_idx,:,:));
    PEAKS{end+1}=data_out.PEAKS(valid_idx,:,:);
    SHIFTS{end+1}=data_out.SHIFTS(valid_idx,:,:);
    SESSION{end+1}=strcat(parts{1},'_',parts{2});
    %SHIFTS = cat(1,SHIFTS,data_out.SHIFTS(valid_idx,:,:));
    SPEED = cat(3,SPEED,data_out.SPEED2);
    REG = cat(1,REG,reg(1));
    FID = cat(1,FID,ones(nnz(valid_idx),1)*iF);
end

%%
nChunks = size(SHIFTS{1},3);
x_vec=zeros(nChunks,16);
for ii=1:16
    x_vec(:,ii)=linspace(100,300,nChunks)+(ii-1)*400;
    
end
x_vec = reshape(x_vec,1,[]);
x_vec=x_vec(1:nChunks);
%%
ops = load_default_opt;
figure

    for iR=1:3
        subplot(1,3,iR)
        
        sel_idx = find(startsWith(REG,regions{iR}));
        trial_idx = 1:10;
        FASTEST=zeros(numel(sel_idx),nChunks);
        SLOWEST = FASTEST;
        ALL=FASTEST;
        DIFF=FASTEST;
        for iRep=1:numel(sel_idx)
            sp=SPEED(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp(trial_idx,:));
            [mi,slowest]=min(sp(trial_idx,:));
            p1=zeros(1,nChunks);
            p2=p1;
            %convert from subselect to 'global idx'
            DIFF(iRep,:)=ma-mi;
            fastest=fastest+trial_idx(1)-1;
            slowest = slowest+trial_idx(1)-1;
            for ii=1:nChunks
                FASTEST(iRep,ii)=nanmean(SHIFTS{sel_idx(iRep)}(:,fastest(ii),ii));
                SLOWEST(iRep,ii)=nanmean(SHIFTS{sel_idx(iRep)}(:,slowest(ii),ii));
                ALL(iRep,ii)=nanmean(nanmean(SHIFTS{sel_idx(iRep)}(:,trial_idx,ii),2),1);
            end
            
        end
        
        %plot(x_vec,nanmean(SLOWEST))
       boundedline(x_vec,nanmean(SLOWEST),nanstd(SLOWEST)/sqrt(size(SLOWEST,1)),'alpha','cmap',[0 0 1])
        hold on
        boundedline(x_vec,nanmean(FASTEST),nanstd(FASTEST)/sqrt(size(SLOWEST,1)),'alpha','cmap',[1 0 0])
        xlabel('window center')
        ylabel('map shift')
        hold on
        %plot(nanmean(ALL))
        title(regions{iR})
        ylim([-4 6])
        yyaxis right
        plot(x_vec,nanmean(DIFF))
        ylim([0 30])
        ylabel('delta speed')
        
        
    end

legend('Slow','Fast','Speed Diff')
legend('Location','SouthWest')
%%
figure
for iR=1:3
        subplot(1,2,1)
        hold on
        sel_idx = find(startsWith(REG,regions{iR}));
        trial_idx = 1:10;
        FASTEST=zeros(numel(sel_idx),nChunks);
        SLOWEST = FASTEST;
        ALL=FASTEST;
        DIFF=FASTEST;
        for iRep=1:numel(sel_idx)
            sp=SPEED(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp(trial_idx,:));
            [mi,slowest]=min(sp(trial_idx,:));
            p1=zeros(1,nChunks);
            p2=p1;
            %convert from subselect to 'global idx'
            DIFF(iRep,:)=ma-mi;
            fastest=fastest+trial_idx(1)-1;
            slowest = slowest+trial_idx(1)-1;
            for ii=1:nChunks
                FASTEST(iRep,ii)=nanmean(SHIFTS{sel_idx(iRep)}(:,fastest(ii),ii));
                SLOWEST(iRep,ii)=nanmean(SHIFTS{sel_idx(iRep)}(:,slowest(ii),ii));
                ALL(iRep,ii)=nanmean(nanmean(SHIFTS{sel_idx(iRep)}(:,trial_idx,ii),2),1);
            end
            
        end
        mm=FASTEST-SLOWEST;
        boundedline(x_vec,nanmean(mm),nanstd(mm)/sqrt(size(mm,1)),'cmap',ops.colors.(regions{iR}),'alpha')
        dd=mm./DIFF;
        subplot(1,2,2)
        hold on
        boundedline(x_vec,nanmean(dd),nanstd(dd)/sqrt(size(dd,1)),'cmap',ops.colors.(regions{iR}),'alpha')
        
        %plot(nanmean(ALL))
end
    subplot(1,2,1)
    legend(regions)
    ylabel('difference in shift')
    xlabel('position')
    subplot(1,2,2)
    ylabel('estimated delay')
    xlabel('Position')
    
%% average over reps and sessions
SESSION_DATA=struct();
SESSION_SPEED=struct();
for iR=1:3
    sel_idx = (startsWith(REG,regions{iR}));
    SHIFTS_THIS = SHIFTS(sel_idx);
    SESSION_THIS = SESSION(sel_idx);
    SPEED_THIS = SPEED(:,:,sel_idx);
    
    [uS,~,idx]=unique(SESSION_THIS);
    SESSION_DATA.(regions{iR})=cell(numel(uS),1);
    for iS=1:numel(uS)
        sel_idx = find(idx==iS);

    FASTEST=zeros(nnz(sel_idx),nChunks);
        SLOWEST = FASTEST;
        Speed_diff = FASTEST;

        for iRep=1:numel(sel_idx)
            sp=SPEED_THIS(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp);
            [mi,slowest]=min(sp);
            Speed_diff(iRep,:)= ma-mi;

            p1=zeros(1,nChunks);
            p2=p1;
            for ii=1:nChunks
                FASTEST(iRep,ii)=nanmean(SHIFTS_THIS{sel_idx(iRep)}(:,fastest(ii),ii));
                SLOWEST(iRep,ii)=nanmean(SHIFTS_THIS{sel_idx(iRep)}(:,slowest(ii),ii));
            end
            
        end
        SESSION_DATA.(regions{iR}){iS}=cat(3,FASTEST,SLOWEST);
        SESSION_SPEED.(regions{iR}){iS}=Speed_diff;
    end
    
end

%%
tracefig =figure();
for iR=1:3
    data_this =SESSION_DATA.(regions{iR});
    av_fast = [];
    av_slow = [];
    av_speed = [];
    for iS=1:size(data_this,1)
        if size(data_this{iS},1)<2
            continue
        end
        av_fast=cat(1,av_fast,squeeze(nanmean(data_this{iS}(:,:,1))));
        av_slow=cat(1,av_slow,squeeze(nanmean(data_this{iS}(:,:,2))));
        av_speed = cat(1,av_speed,squeeze(nanmean(SESSION_SPEED.(regions{iR}){iS})));
    end
    figure(tracefig);
        subplot(2,3,iR)

    boundedline(x_vec,nanmean(av_fast),nanstd(av_fast)/sqrt(size(av_fast,1)),'alpha','cmap',[1 0 0])
    boundedline(x_vec,nanmean(av_slow),nanstd(av_slow)/sqrt(size(av_slow,1)),'alpha','cmap',[0 0 1])
        ylabel('map shift')
    
    ylim([-2 2])
        yyaxis right
        boundedline(x_vec,nanmean(av_speed),nanstd(av_speed)/sqrt(size(av_speed,1)),'alpha','cmap',[.5 .5 .5])
    xlabel('Position')
    title(regions{iR})
    subplot(2,3,iR+3)
    plot(x_vec,av_fast'-av_slow','b')
    
        ylim([-10 10])
        tmp=av_fast-av_slow;
        tmps = nanmean(tmp,2);
        tmpd=nanmean(av_speed,2);
    DELAYS{iR}=tmps./tmpd;
end
        figure;plotSpread(DELAYS)
        xticklabels(regions)