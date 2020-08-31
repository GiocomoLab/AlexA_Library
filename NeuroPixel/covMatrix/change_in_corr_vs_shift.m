matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_reg_200cmChunk/*.mat');
pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');

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


%%
% matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_reg_orig/*.mat');
% STAB = [];
% CORR_VEC=[];
% SHIFT_VEC=[];
% SID = [];
% 
% for iF=1:numel(matfiles)
% data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
% 
% stab = nanmean(nanmean(data_out.corrMat(:,1:6,1:6),2),3);
% nC=numel(stab);
% nT = size(data_out.corrMat,2);
% nChunks = size(data_out.PEAKS,3);
% 
% corr_vec =zeros(nC,nT*nChunks);
% shift_vec=corr_vec;
% for iC=1:nC
%     tmp = reshape(squeeze(data_out.PEAKS(iC,:,:))',1,[]);
%     corr_vec(iC,:)=tmp;
%     tmp = reshape(squeeze(data_out.SHIFTS(iC,:,:))',1,[]);
%     shift_vec(iC,:)=tmp;
% 
% end
% STAB=cat(1,STAB,stab);
% CORR_VEC = cat(1,CORR_VEC,corr_vec);
% SHIFT_VEC = cat(1,SHIFT_VEC,shift_vec);
% 
% SID = cat(1,SID,iF*ones(nC,1));
% 
% end
%%

%%
cluster_idx = CLUSTER_GROUP<3;
pval_t = 0.05;

for iJ=1:3:2
pval=zeros(16,size(CORR_VEC,2));
x1=pval;
for iT=1:16
delta_c = CORR_VEC(cluster_idx,(iT-1)*nChunks+1)-CORR_VEC(cluster_idx,(iT-1)*nChunks+(nChunks-iJ));
for ii=1:size(CORR_VEC,2)
    mdl = fitlm(SHIFT_VEC(cluster_idx,ii),delta_c);

    pp=anova(mdl);
    pval(iT,ii) = pp{1,5};
    x1(iT,ii)=mdl.Coefficients{2,1};

end
end

figure
subplot(1,3,1)
imagesc(pval,[0 0.05])
hold on
for ii=1:16
xline(ii*nChunks,'w')
end
subplot(1,3,2)
imagesc(pval<pval_t)
hold on
for ii=1:16
xline(ii*nChunks,'w')
end

subplot(1,3,3)
tmp = x1;
tmp(pval>pval_t)=0;
imagesc(tmp,[-.01 .01])
hold on
for ii=1:16
xline(ii*nChunks,'w')
end
end
%%
cluster_idx = CLUSTER_GROUP<3;

CV = CORR_VEC(cluster_idx,:);
SV = SHIFT_VEC(cluster_idx,:);

pval=zeros(16,size(CORR_VEC,2));
x1=pval;
for iT=1:16
delta_c = CV(:,(iT-1)*nChunks+1)-CV(:,(iT-1)*nChunks+nChunks);
for ii=1:size(CV,2)
    mdl = fitlm(SV(:,ii),delta_c);

    pp=anova(mdl);
    pval(iT,ii) = pp{1,5};
    x1(iT,ii)=mdl.Coefficients{2,1};

end
end
%%
pval_t = 0.05/size(CV,2);
figure
subplot(2,1,1)
trial = 6
imagesc(pval(trial,:)<pval_t)
hold on
for ii=1:16
xline(ii*nChunks,'w')
end
title('significance')
% subplot(1,3,2)
% imagesc(x1,[-.01 .01])
% hold on
% for ii=1:16
% xline(ii*11,'w')
% end
set(gca,'XTick',[nChunks/2:nChunks:size(CORR_VEC,2)],'XTickLabel',[1:16])


subplot(2,1,2)
tmp = x1;
tmp(pval>pval_t)=0;
imagesc(tmp(trial,:),[-.01 .01])
hold on
for ii=1:16
xline(ii*nChunks,'w')
end
set(gca,'XTick',[nChunks/2:nChunks:size(CORR_VEC,2)],'XTickLabel',[1:16])
title('slope in sig. models')
%%
x_vec=zeros(nChunks,16);
for ii=1:16
x_vec(:,ii)=linspace(100,300,nChunks)+(ii-1)*400;

end

take_idx = [];
for ii=0:15
    take_idx=cat(2,take_idx,[1 11]+(11*ii));
end
take_idx=1:16*nChunks;
x_vec = reshape(x_vec,1,[]);
iT=7;
cluster_idx = CLUSTER_GROUP<3;
delta_c = CORR_VEC(:,(iT-1)*nChunks+1)-CORR_VEC(:,(iT-1)*nChunks+nChunks);
delta_c = delta_c(cluster_idx);
CV = CORR_VEC(cluster_idx,:);
SV = SHIFT_VEC(cluster_idx,:);
SC= SHIFT_CLASSIC(cluster_idx,:);
STABLE = SHIFT_CLASSIC(~cluster_idx,:);
%delta_c = nanmean(nanmean(CORR_MAT(:,1:6,7:10),2),3);
%delta_c = delta_c./mean(CORR_VEC(:,1:66),2);
slices = 3;
data2compare = {};
CT=cbrewer('qual', 'Set2', slices);
N=floor(numel(delta_c)/slices);
[a,b]=sort(delta_c);
figure('Renderer','Painters')
for iS=1:slices
idx = (iS-1)*N+1:iS*N;
idx=b(idx);
tmp_corr = mean(CV(idx,take_idx));
tmp_shift = mean(SV(idx,take_idx));
tmp_s = nanmean(SC(idx,:));
data2compare{iS}=SV(idx,take_idx);
subplot(11,1,[1:5])
hold on
boundedline(x_vec(take_idx),tmp_corr,std(CV(idx,take_idx))/sqrt(nnz(idx)),'cmap',CT(iS,:),'alpha')

%plot(x_vec,tmp_corr,'Color',CT(iS+2,:))
subplot(11,1,[6:10])
hold on
boundedline(x_vec(take_idx),tmp_shift,std(SV(idx,take_idx))/sqrt(nnz(idx)),'cmap',CT(iS,:),'alpha')
%plot(x_vec(floor(nChunks/2):nChunks:end),tmp_s*-1,'Color',CT(iS,:),'LineWidth',2)
%lot(x_vec,tmp_shift,'Color',CT(iS+2,:))
end
%plot(x_vec(floor(nChunks/2):nChunks:end),mean(STABLE)*-1,'k--','LineWidth',2)

subplot(11,1,[1:5])
for ii=1:4
xline((5+ii)*400+1)
end
title('XCorrPeak')
legend({'small','large'})
axis tight
box off
subplot(11,1,[6:10])
for ii=1:4
xline((5+ii)*400+1)
end
legend({'small','large'})
title('Shift')
axis tight
box off
pvals = zeros(1,numel(take_idx));
for ii=1:numel(take_idx)
    pvals(ii)=ranksum(data2compare{1}(:,ii),data2compare{2}(:,ii),'tail','left');
end
subplot(11,1,11)
imagesc(pvals<0.025)
box off
%saveas(gcf,'/Users/attialex/Dropbox/temporary_images/fig6_shiftvsgain.pdf')
%% looking at those that remain stable in trial 1
x_vec=zeros(nChunks,16);
for ii=1:16
x_vec(:,ii)=linspace(100,300,nChunks)+(ii-1)*400;

end

take_idx=1:16*nChunks;
x_vec = reshape(x_vec,1,[]);
cluster_idx = CLUSTER_GROUP<3;
% delta_c_1 = CORR_VEC(:,(iT-1)*nChunks+1)-CORR_VEC(:,(iT-1)*nChunks+nChunks);
% delta_c_1 = delta_c_1(cluster_idx);
gain_1_stab = nanmean(nanmean(CORR_MAT(:,1:6,7),2),3);
ranking = 1:numel(gain_1_stab);
[~,tmpsid]=sort(gain_1_stab);
ranking(tmpsid)=ranking;
ranking = ranking/max(ranking);
iT=8;
delta_c = CORR_VEC(:,(iT-1)*nChunks+1)-CORR_VEC(:,(iT-1)*nChunks+10);

selection_idx = cluster_idx & ranking'>0.5;
delta_c = delta_c(selection_idx);

CV = CORR_VEC(selection_idx,:);
SV = SHIFT_VEC(selection_idx,:);
SC= SHIFT_CLASSIC(selection_idx,:);
%delta_c = nanmean(nanmean(CORR_MAT(:,1:6,7:10),2),3);
%delta_c = delta_c./mean(CORR_VEC(:,1:66),2);
slices = 2;
data2compare = {};
CT=cbrewer('qual', 'Set1', slices);
N=floor(numel(delta_c)/slices);
[a,b]=sort(delta_c);
figure('Renderer','Painters')
for iS=1:slices
idx = (iS-1)*N+1:iS*N;
idx=b(idx);
tmp_corr = mean(CV(idx,take_idx));
tmp_shift = mean(SV(idx,take_idx));
tmp_s = nanmean(SC(idx,:));
data2compare{iS}=SV(idx,take_idx);
subplot(11,1,[1:5])
hold on
boundedline(x_vec(take_idx),tmp_corr,std(CV(idx,take_idx))/sqrt(nnz(idx)),'cmap',CT(iS,:),'alpha')

%plot(x_vec,tmp_corr,'Color',CT(iS+2,:))
subplot(11,1,[6:10])
hold on
boundedline(x_vec(take_idx),tmp_shift,std(SV(idx,take_idx))/sqrt(nnz(idx)),'cmap',CT(iS,:),'alpha')
%plot(x_vec(floor(nChunks/2):nChunks:end),tmp_s*-1,'Color',CT(iS+2,:),'LineWidth',2)
%lot(x_vec,tmp_shift,'Color',CT(iS+2,:))
end

subplot(11,1,[1:5])
for ii=1:4
xline((5+ii)*400+1)
end
title('XCorrPeak')
legend({'small','large'})
axis tight
box off
subplot(11,1,[6:10])
for ii=1:4
xline((5+ii)*400+1)
end
legend({'small','large'})
title('Shift')
axis tight
box off
pvals = zeros(1,numel(take_idx));
for ii=1:numel(take_idx)
    pvals(ii)=ranksum(data2compare{1}(:,ii),data2compare{2}(:,ii),'tail','left');
end
subplot(11,1,11)
imagesc(pvals<0.025)
box off
%saveas(gcf,'/Users/attialex/Dropbox/temporary_images/fig6_shiftvsgain.pdf')


%%
idx = diag(true(1,15),1);
d1=CORR_MAT(:,idx);
s1=SHIFT_MAT(:,idx);
delta_c = d1(:,5);
slices = 2;
CT=cbrewer('seq', 'Blues', slices+2);
N=floor(numel(delta_c)/slices);
[a,b]=sort(delta_c);
figure
for iS=1:slices
idx = (iS-1)*N+1:iS*N;
idx=b(idx);
tmp_corr = mean(d1(idx,:));
tmp_shift = mean(s1(idx,:));
subplot(1,2,1)
hold on
plot(tmp_corr,'Color',CT(iS+2,:))
subplot(1,2,2)
hold on
plot(tmp_shift,'Color',CT(iS+2,:))
end
