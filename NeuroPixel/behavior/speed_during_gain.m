files = dir('/Volumes/T7/attialex/tbtxcorr_MEC_0.8/np*.mat');
%%

pca_data=load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/rep_clusters.mat');


%%
CLU_ID =[];
CM = [];
CM_run = [];
SPEED_MAT = [];
TPT=[];
CM_run=[];
for iF=1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    [~,sn]=fileparts(files(iF).name);
    rep_name = strcat(sn(1:end-1),'rep',sn(end));
    idx = find(strcmp(pca_data.rep_id,rep_name));
    if isempty(idx)
    continue
end
    clu_this = pca_data.rep_cluster(idx);
    CLU_ID = cat(1,CLU_ID,clu_this);
    
    reg_idx = startsWith(data.region,'MEC')';
    stab_idx = nanmean(nanmean(data.corrMat(:,1:6,1:6),2),3)>0.5;
    CM = cat(3,CM,squeeze(nanmean(data.corrMat(reg_idx & stab_idx,:,:))));
    mm = data.trial_speed;
    mm(7:10,:,:)=mm(7:10,:,:)/0.8;
    mm=fillmissing(mm,'linear',2);
    mm=mm/prctile(mm(:),99);
    tmp = squareform(1-pdist(mm,'correlation'));
    CM_run = cat(3,CM_run,tmp);
    SPEED_MAT = cat(3,SPEED_MAT,mm);
    time_per_trial = zeros(1,16);
    for iT=1:16
        tmp_t = min(data.trials)-1+iT;
        time_per_trial(iT)=nnz(data.trials==tmp_t);
    end
    TPT = cat(1,TPT,time_per_trial);
end
%% bl_speed
bl_speed = nanmean(SPEED_MAT(1:6,:,:));
gain_speed = nanmean(SPEED_MAT(7:10,:,:),1);
speed_stable = nanmean(bl_speed(:,:,CLU_ID==1),3);
speed_other = nanmean(bl_speed(:,:,CLU_ID>1),3);
figure
subplot(1,2,1)
plot(speed_stable)
hold on
plot(speed_other)
subplot(1,2,2)
tmp=gain_speed./bl_speed;
boundedline(1:2:400,nanmean(tmp(:,:,CLU_ID==1),3),nanstd(tmp(:,:,CLU_ID==1),[],3)/sqrt(nnz(CLU_ID==1)),'alpha')
hold on
boundedline(1:2:400,nanmean(tmp(:,:,CLU_ID>1),3),nanstd(tmp(:,:,CLU_ID>1),[],3)/sqrt(nnz(CLU_ID>1)),'alpha','cmap',[1 0 0])
legend({'Stable Cluster','remap cluster'})
pvals = zeros(1,200);
for iB = 1:200
    pvals(iB)=ranksum(squeeze(tmp(:,iB,CLU_ID==1)),squeeze(tmp(:,iB,CLU_ID>1)));
end
imagesc(1:2:400,0,pvals<0.05)
ylim([.25 1.6])
%% time per trial
tptN=TPT./mean(TPT(:,1:6),2);
figure
boundedline(1:16,mean(tptN(CLU_ID==1,:)),std(tptN(CLU_ID==1,:))/sqrt(nnz(CLU_ID==1)),'alpha')
boundedline(1:16,mean(tptN(CLU_ID>1,:)),std(tptN(CLU_ID>1,:))/sqrt(nnz(CLU_ID>1)),'cmap',[1 0 0],'alpha')

for iT=1:16
    p = ranksum(tptN(CLU_ID>1,iT),tptN(CLU_ID==1,iT));
    text(iT,2,sprintf('%.3e',p))
end
%%
figure
subplot(2,2,1)
imagesc(squeeze(nanmean(CM(:,:,CLU_ID==1),3)),[0 0.7])
axis image
title('spatial sim stable')
subplot(2,2,2)
imagesc(squeeze(nanmean(CM(:,:,CLU_ID>1),3)),[0 0.7])
axis image
title('spatial sim remap')
subplot(2,2,3)
imagesc(squeeze(nanmean(CM_run(:,:,CLU_ID==1),3)),[0 0.7])
axis image
title('running sim stable')


subplot(2,2,4)
imagesc(squeeze(nanmean(CM_run(:,:,CLU_ID>1),3)),[0 0.7])
axis image
title('runing sim remap')