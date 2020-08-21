matfiles = dir('/Users/attialex/slidingWindow_baseline_100cm/*.mat');

regions = {'MEC','VISp','RS'};
PEAKS = {};
SHIFTS = {};
SPEED = [];
REG={};
CLU=[];
FID = [];
for iF=1:numel(matfiles)
    
    [~,sn]=fileparts(matfiles(iF).name);
    parts = strsplit(sn,'_');
    if str2double(parts{end})>5
        continue
    end
    
    
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
    reg=reg(valid_idx);
    %PEAKS = cat(1,PEAKS,data_out.PEAKS(valid_idx,:,:));
    PEAKS{end+1}=data_out.PEAKS(valid_idx,:,:);
    SHIFTS{end+1}=data_out.SHIFTS(valid_idx,:,:);
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
    
%%
