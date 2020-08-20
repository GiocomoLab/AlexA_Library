pca_data = load('G:\My Drive\rep_clusters.mat');
matfiles = dir('F:\attialex\slidingWindow_08\*.mat');

regions = {'MEC','VISp','RS'};
PEAKS = {};
SHIFTS = {};
SPEED = [];
REG={};
CLU=[];
FID = [];
for iF=1:numel(matfiles)
    
    [~,sn]=fileparts(matfiles(iF).name);
    tag = sn(end);
    sn_rep = strcat(sn(1:end-1),'rep',tag);
    if ~ismember(sn_rep, pca_data.rep_id)
        continue
    else
        idx = find(startsWith(pca_data.rep_id,sn_rep));
        stab_cluster = pca_data.rep_cluster(idx);
        reg_cluster = pca_data.brain_region(idx);
    end
    
    
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    nC=size(data_out.corrMat,1);
    stab = nanmean(nanmean(data_out.corrMat(:,1:6,1:6),2),3);
    reg = data_out.region;
    idx = startsWith(reg,'VISpm');
    reg(idx)={'VISm'};
    idx = startsWith(reg,'RSPag');
    reg(idx)={'aglRS'};
    idx = startsWith(reg,'RS');
    reg(idx)={'RSC'};
    
    valid_idx = stab'>=0.5 & startsWith(reg,regions);
    if nnz(valid_idx)<5
        disp
        continue
    end
    
    %PEAKS = cat(1,PEAKS,data_out.PEAKS(valid_idx,:,:));
    PEAKS{end+1}=data_out.PEAKS(valid_idx,:,:);
    SHIFTS{end+1}=data_out.SHIFTS(valid_idx,:,:);
    %SHIFTS = cat(1,SHIFTS,data_out.SHIFTS(valid_idx,:,:));
    SPEED = cat(3,SPEED,data_out.SPEED2);
    CLU = cat(1,CLU,stab_cluster);
    REG = cat(1,REG,reg_cluster);
    FID = cat(1,FID,ones(nnz(valid_idx),1)*iF);
end

%%
nChunks = 21;
x_vec=zeros(nChunks,16);
for ii=1:16
    x_vec(:,ii)=linspace(100,300,nChunks)+(ii-1)*400;
    
end
x_vec = reshape(x_vec,1,[]);
x_vec=x_vec(1:nChunks);
%%
figure
trial_block={[1:6],[7:10],[11:16]};
for iTB=1:3
    for iR=1:3
        subplot(3,3,iR+(iTB-1)*3)
        
        sel_idx = find(startsWith(REG,regions{iR}) & CLU==1);
        trial_idx = trial_block{iTB};
        FASTEST=zeros(numel(sel_idx),21);
        SLOWEST = FASTEST;
        ALL=FASTEST;
        DIFF=FASTEST;
        for iRep=1:numel(sel_idx)
            sp=SPEED(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp(trial_idx,:));
            [mi,slowest]=min(sp(trial_idx,:));
            p1=zeros(1,21);
            p2=p1;
            %convert from subselect to 'global idx'
            DIFF(iRep,:)=ma-mi;
            fastest=fastest+trial_idx(1)-1;
            slowest = slowest+trial_idx(1)-1;
            for ii=1:21
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
end
legend('Slow','Fast','Speed Diff')
legend('Location','SouthWest')
%%
figure
trial_block={[1:6],[7:10],[11:16]};
for iR=1:3
    for iTB=1:3
        subplot(1,3,iR)
        hold on
        sel_idx = find(startsWith(REG,regions{iR}) & CLU==1);
        trial_idx = trial_block{iTB};
        FASTEST=zeros(numel(sel_idx),21);
        SLOWEST = FASTEST;
        ALL=FASTEST;
        DIFF=FASTEST;
        for iRep=1:numel(sel_idx)
            sp=SPEED(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp(trial_idx,:));
            [mi,slowest]=min(sp(trial_idx,:));
            p1=zeros(1,21);
            p2=p1;
            %convert from subselect to 'global idx'
            DIFF(iRep,:)=ma-mi;
            fastest=fastest+trial_idx(1)-1;
            slowest = slowest+trial_idx(1)-1;
            for ii=1:21
                FASTEST(iRep,ii)=nanmean(SHIFTS{sel_idx(iRep)}(:,fastest(ii),ii));
                SLOWEST(iRep,ii)=nanmean(SHIFTS{sel_idx(iRep)}(:,slowest(ii),ii));
                ALL(iRep,ii)=nanmean(nanmean(SHIFTS{sel_idx(iRep)}(:,trial_idx,ii),2),1);
            end
            
        end
        
        plot(x_vec,nanmean(FASTEST-SLOWEST))
        %plot(nanmean(ALL))
    end
    title(regions{iR})
    legend('BL PRE','GAIN','BL POST')
    ylim([-2 5])
end
%%
%%
figure
trial_block={[1:6],[7:10],[11:16]};
for iTB=1:3
    for iR=1:3
        subplot(3,3,iR+(iTB-1)*3)
        
        sel_idx = find(startsWith(REG,regions{iR}) & CLU==1);
        trial_idx = trial_block{iTB};
        FASTEST=zeros(numel(sel_idx),21);
        SLOWEST = FASTEST;
        ALL=FASTEST;
        DIF=FASTEST;
        for iRep=1:numel(sel_idx)
            sp=SPEED(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp(trial_idx,:));
            [mi,slowest]=min(sp(trial_idx,:));
            p1=zeros(1,21);
            p2=p1;
            %convert from subselect to 'global idx'
            DIFF(iRep,:)=ma-mi;
            fastest=fastest+trial_idx(1)-1;
            slowest = slowest+trial_idx(1)-1;
            for ii=1:21
                FASTEST(iRep,ii)=nanmean(PEAKS{sel_idx(iRep)}(:,fastest(ii),ii));
                SLOWEST(iRep,ii)=nanmean(PEAKS{sel_idx(iRep)}(:,slowest(ii),ii));
                ALL(iRep,ii)=nanmean(nanmean(PEAKS{sel_idx(iRep)}(:,trial_idx,ii),2),1);
            end
            
        end
        
        %plot(x_vec,nanmean(SLOWEST))
        boundedline(x_vec,nanmean(SLOWEST),nanstd(SLOWEST)/sqrt(size(SLOWEST,1)),'alpha')
        hold on
        boundedline(x_vec,nanmean(FASTEST),nanstd(FASTEST)/sqrt(size(SLOWEST,1)),'alpha')

        plot(x_vec,nanmean(FASTEST))
        hold on
        %plot(nanmean(ALL))
        title(regions{iR})
        
        
        
        
    end
end
legend('Slow','Fast','Speed Diff')
legend('Location','SouthWest')
%% fastest-slowest averaged across cells only

figure
trial_block={[1:6],[7:10],[11:16]};
for iR=1:3
    for iTB=1:3
        subplot(1,3,iR)
        hold on
        sel_idx = find(startsWith(REG,regions{iR}) & CLU==1);
        trial_idx = trial_block{iTB};
        DIFF=[];
        for iRep=1:numel(sel_idx)
            sp=SPEED(:,:,sel_idx(iRep));
            [ma,fastest]=max(sp(trial_idx,:));
            [mi,slowest]=min(sp(trial_idx,:));
            p1=zeros(1,21);
            p2=p1;
            %convert from subselect to 'global idx'
            
            fastest=fastest+trial_idx(1)-1;
            slowest = slowest+trial_idx(1)-1;
            tmpf=zeros(size(SHIFTS{sel_idx(iRep)},1),21);
            tmps=tmpf;
            for ii=1:21
                tmpf(:,ii) = SHIFTS{sel_idx(iRep)}(:,fastest(ii),ii);
                tmps(:,ii) = SHIFTS{sel_idx(iRep)}(:,slowest(ii),ii);
                
                
            end
            DIFF=cat(1,DIFF,tmpf-tmps);
        end
        
        plot(x_vec,nanmean(DIFF))
        hold on
        title(regions{iR})
            ylim([-2 5])

        
        
        
    end
end
title(regions{iR})
    legend('BL PRE','GAIN','BL POST')
    ylim([-2 5])

