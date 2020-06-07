ops.towerbins = -41:2:41;
ops.towerbincent = ops.towerbins(1:end-1)*.5 + ops.towerbins(2:end)*.5;
ops.xbincent = 1:2:400;
regions = {'MEC','VISp','RS'};
files = dir('/Volumes/Samsung_T5/attialex/pca_classifier_all/*.mat');
files= dir('/Volumes/Samsung_T5/attialex/circular_classifier_3/*.mat');
%%

AVG_ERR = struct();
AVG_ERR_MLD = struct();
AVG_ERR_TOWER = struct();
AVG_ERR_TOWER_MLD = struct();
for iR=1:numel(regions)
    AVG_ERR.(regions{iR})=[];
    AVG_ERR_MLD.(regions{iR})=[];
    AVG_ERR_TOWER.(regions{iR})=[];
    
    AVG_ERR_TOWER_MLD.(regions{iR})=[];
end

for iF=1:numel(files)
    data_out = load(fullfile('/Volumes/Samsung_T5/attialex/circular_classifier_3',files(iF).name));
    
    
    tower_relative = data_out.true_pos;
    tower_relative = mod(tower_relative+40,80) - 40 +1; %now relative to each of the
    tower_relative_bin = discretize(tower_relative,ops.towerbins);
    
    [~,~,pos_bin] = histcounts(data_out.true_pos,0:2:400);
    pos_bin(pos_bin==0)=1;
    
    true_pos = ops.xbincent(pos_bin);
    
    avg_cd = zeros(1,200);
    avg_mld = avg_cd;
    
%     err = data_out.error_cd;
%     correction_idx = abs(err)>pi;
%     err(correction_idx) = err(correction_idx)-pi*sign(err(correction_idx));
%     err_cd = err;
%     err_cd(abs(err_cd)>pi/4)=nan;
%     err_cd = sin(err_cd);
    err= data_out.true_theta-data_out.predicted_theta;
    score = cos(err);
    adj_score = (acos(score)/-pi +.5)*2;
    %err_cd = adj_score;
    err_cd = (1-(acos(adj_score)/-pi+.5)*2)*100.*sign(err);
    err_cd(abs(err_cd)>50)=nan;
    err = ops.xbincent(data_out.pos_bin+1)-ops.xbincent(data_out.predicted_bin);
    correction_idx = abs(err)>400/2;
    err(correction_idx) = err(correction_idx)-400*sign(err(correction_idx));
    err_mld = err;
    err_mld(abs(err_mld)>50)=nan;
    for ii=1:200
        idx = pos_bin==ii;
        avg_cd(ii)=nanmean(err_cd(idx));
        avg_mld(ii)=nanmean(err_mld(idx));
    end
    relative_tower_error_circular=zeros(1,numel(ops.towerbincent));
    relative_tower_error_mld=zeros(1,numel(ops.towerbincent));

    for ij = 1:numel(ops.towerbincent)
        idx = tower_relative_bin==ij & pos_bin<180 & pos_bin>20;
        relative_tower_error_circular(ij)=nanmean(err_cd(idx));
        relative_tower_error_mld(ij)=nanmean(err_mld(idx));

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
    AVG_ERR.(winner) = cat(1,AVG_ERR.(winner),avg_cd);
    AVG_ERR_MLD.(winner) = cat(1,AVG_ERR_MLD.(winner),avg_mld);
    AVG_ERR_TOWER.(winner)=cat(1,AVG_ERR_TOWER.(winner),relative_tower_error_circular);
    AVG_ERR_TOWER_MLD.(winner)=cat(1,AVG_ERR_TOWER_MLD.(winner),relative_tower_error_mld);

end
%%
figure;
hold on;
cmap=cbrewer('qual','Set2',3,'pchip');
for iR=1:3
    
    subplot(1,2,1)
    hold on
    nn = nanmean(AVG_ERR.(regions{iR}));
    ee = nanstd(AVG_ERR.(regions{iR}))/sqrt(size(AVG_ERR.(regions{iR}),1));
    
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    title('Circular')
    subplot(1,2,2)
    hold on
    nn=nanmean(AVG_ERR_MLD.(regions{iR}));
    ee=nanstd(AVG_ERR_MLD.(regions{iR}))/sqrt(size(AVG_ERR.(regions{iR}),1));
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    title('MLD')
    %plot(nn)
    
end
subplot(1,2,1)
legend(regions)
xlabel('Position')
ylabel('error')
subplot(1,2,2)
legend(regions)
xlabel('Position')
ylabel('error')
%%

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
    title('Circular')

       subplot(1,2,2)
    hold on
    nn=nanmean(AVG_ERR_TOWER_MLD.(regions{iR}));
    ee=nanstd(AVG_ERR_TOWER_MLD.(regions{iR}))/sqrt(size(AVG_ERR.(regions{iR}),1));
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',cmap(iR,:))
    title('MLD')
    
end
subplot(1,2,1)
legend(regions)
xlabel('Distance to tower')
ylabel('error')

subplot(1,2,2)
legend(regions)
xlabel('Distance to tower')
ylabel('error')

%% TOWER 2
figure;
hold on;
cmap=cbrewer('qual','Set2',3,'pchip');
for iR=1:3
    
    subplot(1,3,iR)
    hold on
    tmp = AVG_ERR_MLD.(regions{iR});
    tmp = cat(3,tmp(:,20:60),tmp(:,60:100),tmp(:,100:140),tmp(:,140:180));
    tmp = nanmean(tmp,3);
    nn = nanmean(tmp);
    ee = nanstd(tmp)/sqrt(size(tmp,1));
    
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',[0 0 1])
    
    hold on
    tmp = AVG_ERR.(regions{iR});
    tmp = cat(3,tmp(:,20:60),tmp(:,60:100),tmp(:,100:140),tmp(:,140:180));
    %tmp = cat(3,tmp(:,60:100),tmp(:,140:180));
    tmp = nanmean(tmp,3);
    nn = nanmean(tmp);
    ee = nanstd(tmp)/sqrt(size(tmp,1));
    
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',[1 0 0])
    title(regions{iR})
    legend('FiringRate','Circular')
    xlabel('Distance from Tower')
    
end
