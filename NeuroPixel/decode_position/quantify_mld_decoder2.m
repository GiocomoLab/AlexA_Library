%% evaluate classifier
ops = load_default_opt;

ops.towerbins = -41:2:41;
ops.towerbincent = ops.towerbins(1:end-1)*.5 + ops.towerbins(2:end)*.5;
ops.xbin = 2;
ops.xbinedges = 0:ops.xbin:400;
ops.xbincent = ops.xbinedges(1:end-1)+ops.xbin/2;
ops.nBins = numel(ops.xbincent);
ops.edges = ops.xbinedges;
ops.filter = gausswin(11);
ops.trials=5:20;
ops.num_pcs_decoding = 10;
regions = {'MEC','VISp','RS'};

regions = {'MEC','VISp','RS'};
files = dir('/Volumes/Samsung_T5/attialex/mld_classifier2/*.mat');
AVG_ERR = struct();
AVG_ABSERR =struct();
AVG_ABSERR_TOWER = struct();
AVG_ERR_TOWER = struct();
AVG_CM = struct();
for iR=1:numel(regions)
    AVG_ERR.(regions{iR})=[];
    AVG_ABSERR.(regions{iR})=[];
    AVG_ABSERR_TOWER.(regions{iR})=[];
    AVG_ERR_TOWER.(regions{iR})=[];
    AVG_CM.(regions{iR})=[];
end

for iF=1:numel(files)
    data_out = load(fullfile(files(iF).folder,files(iF).name));
    
    
    tower_relative = repmat(ops.xbincent,1,16);
    tower_relative = mod(tower_relative+40,80) - 40 +1; %now relative to each of the
    tower_relative_bin = discretize(tower_relative,ops.towerbins);
    pos_bin = repmat(1:200,1,16);
    [N,Xedges,Yedges] = histcounts2(ops.xbincent(pos_bin),data_out.pred_pos,ops.xbinedges,ops.xbinedges,'Normalization','probability');
    true_pos = ops.xbincent(pos_bin);
    err = data_out.pred_pos-true_pos;
    correction_idx = abs(err)>400/2;
    err(correction_idx) = err(correction_idx)-400*sign(err(correction_idx));
    err_this = err;
    err_this(abs(err_this)>50)=nan;
    
    abs_err = abs(err_this);
   
    avg_err = zeros(1,200);
    avg_abs_err = avg_err;
    for ii=1:200
        idx = pos_bin==ii;
        avg_err(ii)=nanmean(err_this(idx));
        avg_abs_err(ii)=nanmean(abs_err(idx));
        
    end
    relative_tower_error=zeros(1,numel(ops.towerbincent));
    relative_tower_abs_error = relative_tower_error;
    for ij = 1:numel(ops.towerbincent)
        idx = tower_relative_bin==ij & pos_bin<180 & pos_bin>20;
        relative_tower_error(ij)=nanmean(err_this(idx));
        relative_tower_abs_error(ij)=nanmean(abs_err(idx));
    end
    
    counts = nan(size(regions));
    for iR=1:3
        counts(iR)=nnz(startsWith(data_out.region,regions{iR}));
    end
    [a,ii]=max(counts);
    winner = regions{ii};
    if startsWith(winner,'RS')
        winner = 'RS';
    end
    AVG_ERR.(winner) = cat(1,AVG_ERR.(winner),avg_err);
    AVG_ABSERR.(winner)= cat(1,AVG_ABSERR.(winner),avg_abs_err);
    AVG_ERR_TOWER.(winner)=cat(1,AVG_ERR_TOWER.(winner),relative_tower_error);
    AVG_ABSERR_TOWER.(winner) = cat(1,AVG_ABSERR_TOWER.(winner),relative_tower_abs_error);
    AVG_CM.(winner) = cat(3,AVG_CM.(winner),N);
end
%%
figure;
hold on;
cmap=cbrewer('qual','Set2',3,'pchip');
for iR=1:3
    
    subplot(2,1,1)
    hold on
    nn = nanmean(AVG_ERR.(regions{iR}));
    ee = nanstd(AVG_ERR.(regions{iR}))/sqrt(size(AVG_ERR.(regions{iR}),1));
    
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    
       
    subplot(2,1,2)
    hold on
    nn = nanmean(AVG_ABSERR.(regions{iR}));
    ee = nanstd(AVG_ABSERR.(regions{iR}))/sqrt(size(AVG_ABSERR.(regions{iR}),1));
    
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    
    
end
subplot(2,1,1)
legend(regions)
xlabel('Position')
ylabel('error')

subplot(2,1,2)
legend(regions)
xlabel('Position')
ylabel('|error|')


%% TOWER 1
figure;
hold on;
cmap=cbrewer('qual','Set2',3,'pchip');
for iR=1:3
    
    subplot(1,2,1)
    hold on
    nn = nanmean(AVG_ERR_TOWER.(regions{iR}));
    ee = nanstd(AVG_ERR_TOWER.(regions{iR}))/sqrt(size(AVG_ERR.(regions{iR}),1));
    
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    
    subplot(1,2,2)
    hold on
    nn=nanmean(AVG_ABSERR_TOWER.(regions{iR}));
    ee=nanstd(AVG_ABSERR_TOWER.(regions{iR}))/sqrt(size(AVG_ABSERR_TOWER.(regions{iR}),1));
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    %plot(nn)
    
end
subplot(1,2,1)
legend(regions)
xlabel('Position')
ylabel('error')
subplot(1,2,2)
legend(regions)
xlabel('Position')
ylabel('|error|')

%% TOWER 2
figure;
hold on;
cmap=cbrewer('qual','Set2',3,'pchip');
for iR=1:3
    
    subplot(1,2,1)
    hold on
    tmp = AVG_ERR.(regions{iR});
    tmp = cat(3,tmp(:,20:60),tmp(:,60:100),tmp(:,100:140),tmp(:,140:180));
    tmp = nanmean(tmp,3);
    nn = nanmean(tmp);
    ee = nanstd(tmp)/sqrt(size(tmp,1));
    
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    
    subplot(1,2,2)
    hold on
    tmp = AVG_ABSERR.(regions{iR});
    tmp = cat(3,tmp(:,20:60),tmp(:,60:100),tmp(:,100:140),tmp(:,140:180));
    tmp = nanmean(tmp,3);
    nn = nanmean(tmp);
    ee = nanstd(tmp)/sqrt(size(tmp,1));
    
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    
end
subplot(1,2,1)
legend(regions)
xlabel('Position')
ylabel('error')

%%
figure
for iR=1:numel(regions)
    subplot(1,3,iR)
    imagesc((mean(AVG_CM.(regions{iR}),3)),[0 0.0001])
    axis image
    title(regions{iR})
end
    colormap(cbrewer('seq','BuPu',20))