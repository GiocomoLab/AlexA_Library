%% evaluate classifier
ops.towerbins = -41:2:41;
ops.towerbincent = ops.towerbins(1:end-1)*.5 + ops.towerbins(2:end)*.5;

regions = {'MEC','VISp','RS'};
files = dir('/Volumes/Samsung_T5/attialex/pca_classifier_TreeBagger/*.mat');
AVG_ERR = struct();
AVG_ERR_SVM = struct();
AVG_ERR_TOWER = struct();
for iR=1:numel(regions)
    AVG_ERR.(regions{iR})=[];
    AVG_ERR_SVM.(regions{iR})=[];
    AVG_ERR_TOWER.(regions{iR})=[];
end

for iF=1:numel(files)
    data_out = load(fullfile(files(iF).folder,files(iF).name));
    
    pred_pos = [];
    pred_pos_svm = [];
    for iFold = 1:numel(data_out.pred_pos)
        pred_pos = cat(2,pred_pos,data_out.pred_pos{iFold});
        pred_pos_svm = cat(2,pred_pos_svm,data_out.yhat{iFold});
    end
    tower_relative = data_out.true_pos;
    tower_relative = mod(tower_relative+40,80) - 40 +1; %now relative to each of the
    tower_relative_bin = discretize(tower_relative,ops.towerbins);
    [~,~,pos_bin] = histcounts(data_out.true_pos,0:2:400);
    pos_bin(pos_bin==0)=1;
    true_pos = ops.xbincent(pos_bin);
    avg_err = zeros(1,200);
    avg_svm = avg_err;
    err = pred_pos-true_pos;
    correction_idx = abs(err)>400/2;
    err(correction_idx) = err(correction_idx)-400*sign(err(correction_idx));
    err_this = err;
    err_this(abs(err_this)>20)=nan;
    
    err = pred_pos_svm-true_pos;
    correction_idx = abs(err)>400/2;
    err(correction_idx) = err(correction_idx)-400*sign(err(correction_idx));
    err_svm = err;
    err_svm(abs(err_svm)>20)=nan;
    
    for ii=1:200
        idx = pos_bin==ii;
        avg_err(ii)=nanmean(err_this(idx));
        avg_svm(ii)=nanmean(err_svm(idx));
    end
    relative_tower_error=zeros(1,numel(ops.towerbincent));
    for ij = 1:numel(ops.towerbincent)
        idx = tower_relative_bin==ij & pos_bin<180 & pos_bin>20;
        relative_tower_error(ij)=nanmean(err_this(idx));
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
    AVG_ERR_SVM.(winner) = cat(1,AVG_ERR_SVM.(winner),avg_svm);
    AVG_ERR_TOWER.(winner)=cat(1,AVG_ERR_TOWER.(winner),relative_tower_error);
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
    
    subplot(1,2,2)
    hold on
    nn=nanmean(AVG_ERR_SVM.(regions{iR}));
    ee=nanstd(AVG_ERR_SVM.(regions{iR}))/sqrt(size(AVG_ERR.(regions{iR}),1));
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',cmap(iR,:))
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
    nn=nanmean(AVG_ERR_SVM.(regions{iR}));
    ee=nanstd(AVG_ERR_SVM.(regions{iR}));
    %boundedline(ops.xbincent,nn,ee)
    plot(nn)
    
end
subplot(1,2,1)
legend(regions)
xlabel('Position')
ylabel('error')

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
    tmp = AVG_ERR_SVM.(regions{iR});
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