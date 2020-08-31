%% load 08 data
matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_reg_200cmChunk/*.mat');
pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');
nChunks = 21;
x_vec=zeros(nChunks,16);
for ii=1:16
x_vec(:,ii)=linspace(100,300,nChunks)+(ii-1)*400;

end
x_vec = reshape(x_vec,1,[]);

STAB = [];
CORR_VEC=[];
SHIFT_VEC=[];
SHIFT_CLASSIC = [];
SID = [];
CORR_MAT =[];
CLUSTER_GROUP = [];
for iF=1:numel(matfiles)
data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
tmp_idx = (strcmp(pca_data.MEC_NAMES,matfiles(iF).name(1:end-4)));
cluster_this = pca_data.MEC_CLUSTERS(tmp_idx);
stab = nanmean(nanmean(data_out.corrMat(:,1:6,1:6),2),3);
shift_classic = squeeze(nanmean(data_out.shiftMat(:,1:6,:),2));
reg_idx = startsWith(data_out.region,'MEC');
if isrow(reg_idx)
reg_idx = reg_idx';
end
idx =stab>=.5 & reg_idx;

if nnz(idx)<5
    continue
end
nC=numel(stab);
nT = size(data_out.corrMat,2);
nChunks = size(data_out.PEAKS,3);

corr_vec =zeros(nC,nT*nChunks);
shift_vec=corr_vec;
PEAKS = data_out.PEAKS;
SHIFTS = data_out.SHIFTS;
%valid_idx = ismember(SHIFTS,[-30 30]);
%SHIFTS(valid_idx)=nan;
%PEAKS(valid_idx)=nan;

for iC=1:nC
    tmp = reshape(squeeze(PEAKS(iC,:,:))',1,[]);
    corr_vec(iC,:)=tmp;
    tmp = reshape(squeeze(SHIFTS(iC,:,:))',1,[]);
    shift_vec(iC,:)=tmp;

end

CORR_VEC = cat(1,CORR_VEC,nanmean(corr_vec(idx,:)));
SHIFT_VEC = cat(1,SHIFT_VEC,nanmean(shift_vec(idx,:)));
SHIFT_CLASSIC = cat(1,SHIFT_CLASSIC,nanmean(shift_classic(idx,:)));
CORR_MAT = cat(1,CORR_MAT,nanmean(data_out.corrMat(idx,:,:)));
CLUSTER_GROUP = cat(1,CLUSTER_GROUP,cluster_this);
end

%% load 05 data
matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_05_reg_200cmChunk/*.mat');

CORR_VEC_05=[];
SHIFT_VEC_05=[];
for iF=1:numel(matfiles)
data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
tmp_idx = (strcmp(pca_data.MEC_NAMES,matfiles(iF).name(1:end-4)));
cluster_this = pca_data.MEC_CLUSTERS(tmp_idx);
stab = nanmean(nanmean(data_out.corrMat(:,1:6,1:6),2),3);
reg_idx = startsWith(data_out.region,'MEC');
if isrow(reg_idx)
reg_idx = reg_idx';
end
idx =stab>=.5 & reg_idx;

if nnz(idx)<5
    continue
end
nC=numel(stab);
nT = size(data_out.corrMat,2);
nChunks = size(data_out.PEAKS,3);

corr_vec =zeros(nC,nT*nChunks);
shift_vec=corr_vec;
PEAKS = data_out.PEAKS;
SHIFTS = data_out.SHIFTS;
%valid_idx = ismember(SHIFTS,[-30 30]);
%SHIFTS(valid_idx)=nan;
%PEAKS(valid_idx)=nan;

for iC=1:nC
    tmp = reshape(squeeze(PEAKS(iC,:,:))',1,[]);
    corr_vec(iC,:)=tmp;
    tmp = reshape(squeeze(SHIFTS(iC,:,:))',1,[]);
    shift_vec(iC,:)=tmp;

end

CORR_VEC_05 = cat(1,CORR_VEC_05,nanmean(corr_vec(idx,:)));
SHIFT_VEC_05 = cat(1,SHIFT_VEC_05,nanmean(shift_vec(idx,:)));

end
%% plot different cluster groups
CT=cbrewer('qual', 'Set2', 3);
figure
for iG=1:2
    cluster_idx = CLUSTER_GROUP==iG;
    tmp_corr = mean(CORR_VEC(cluster_idx,:));
tmp_shift = mean(SHIFT_VEC(cluster_idx,:));
subplot(2,1,1)
hold on
boundedline(x_vec,tmp_corr,std(CORR_VEC(cluster_idx,:))/sqrt(nnz(cluster_idx)),'cmap',CT(iG,:),'alpha')

%plot(x_vec,tmp_corr,'Color',CT(iS+2,:))
subplot(2,1,2)
hold on
boundedline(x_vec,tmp_shift,std(SHIFT_VEC(cluster_idx,:))/sqrt(nnz(cluster_idx)),'cmap',CT(iG,:),'alpha')
end
for iS=1:2
    subplot(2,1,iS)
for ii=1:4
xline((5+ii)*400+1)
end
end
 subplot(2,1,2)
 title('shift')
% boundedline(x_vec,mean(SHIFT_VEC_05),std(SHIFT_VEC_05)/sqrt(size(SHIFT_VEC_05,1)),'cmap',[.4 .4 .4])
 subplot(2,1,1)
 title('similarity')
% boundedline(x_vec,mean(CORR_VEC_05),std(CORR_VEC_05)/sqrt(size(SHIFT_VEC_05,1)),'cmap',[.4 .4 .4])
%%
CT=cbrewer('qual', 'Set2', 3);
figure('Color','white','Renderer','Painters')

    cluster_idx = CLUSTER_GROUP<3;
    tmp_corr = mean(CORR_VEC(cluster_idx,:));
tmp_shift = mean(SHIFT_VEC(cluster_idx,:));
subplot(2,1,1)
hold on
boundedline(x_vec,tmp_corr,std(CORR_VEC(cluster_idx,:))/sqrt(nnz(cluster_idx)),'cmap',CT(1,:),'alpha')

%plot(x_vec,tmp_corr,'Color',CT(iS+2,:))
subplot(2,1,2)
hold on
boundedline(x_vec,tmp_shift,std(SHIFT_VEC(cluster_idx,:))/sqrt(nnz(cluster_idx)),'cmap',CT(1,:),'alpha')

    cluster_idx = CLUSTER_GROUP==3;
    tmp_corr = mean(CORR_VEC(cluster_idx,:));
tmp_shift = mean(SHIFT_VEC(cluster_idx,:));
subplot(2,1,1)
hold on
boundedline(x_vec,tmp_corr,std(CORR_VEC(cluster_idx,:))/sqrt(nnz(cluster_idx)),'cmap',CT(2,:),'alpha')

%plot(x_vec,tmp_corr,'Color',CT(iS+2,:))
subplot(2,1,2)
hold on
boundedline(x_vec,tmp_shift,std(SHIFT_VEC(cluster_idx,:))/sqrt(nnz(cluster_idx)),'cmap',CT(2,:),'alpha')

for iS=1:2
    subplot(2,1,iS)
for ii=1:15
xline((ii)*400,'--')
end
end
subplot(2,1,2)
title('shift')
boundedline(x_vec,mean(SHIFT_VEC_05),std(SHIFT_VEC_05)/sqrt(size(SHIFT_VEC_05,1)),'cmap',[.4 .4 .4])
subplot(2,1,1)
title('similarity')
boundedline(x_vec,mean(CORR_VEC_05),std(CORR_VEC_05)/sqrt(size(SHIFT_VEC_05,1)),'cmap',[.4 .4 .4])
set(gca,'XTick',[200:400:32000],'XTickLabel',[1:16])

xlim([0 4000])
ax1=gca;
ax1.XAxis.TickLength = [0 0];
legend({'Stable 0.8','Remap 0.8', '0.5'})
saveas(gcf,'/Users/attialex/Dropbox/temporary_images/fig6_sliding_gain05.pdf')
