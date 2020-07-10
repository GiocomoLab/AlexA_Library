matfiles = dir('/Volumes/Samsung_T5/attialex/autocorr_08/*.mat');
pca_data=load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');
%%
xcorr_path = '/Volumes/Samsung_T5/tbtxcorr_contigFR_30cm_data_corrected2';
XC_BL=[];
XC_G=[];
XC_G2=[];
Cluster_Group=[];
x_vec = -200:2:200;

zero_bin = find(x_vec==0);

PEAKS_BL  = [];
PEAKS_GC = [];
TROUGHS_BL = [];
TROUGHS_GC = [];
for iF=1:numel(matfiles)
    data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    data_xcorr = load(fullfile(xcorr_path,matfiles(iF).name));
    stab=nanmean(nanmean(data_xcorr.corrMat(:,1:6,1:6),2),3);
    stab_idx=stab>=0.5;
    if nnz(stab_idx)<5
        continue
    end
    xc_bl = nanmean(data.xc_bl(stab_idx,:));
    xc_g = nanmean(data.xc_gain(stab_idx,:));
    %xc_g2 = nanmean(data.xc_gain_corrected(stab_idx,:));
    [~,sn]=fileparts(matfiles(iF).name);
    MEC_idx = find(startsWith(pca_data.MEC_NAMES,sn));
    if ~isempty(MEC_idx)
        cg=pca_data.MEC_CLUSTERS(MEC_idx);
    else
        cg= nan;
    end
    
    nC=numel(stab);
    peaks = nan(nC,1);
    peaks_gc = peaks;
    trough=peaks;
    for iC=1:nC
        tmp = data.xc_bl(iC,:);
        tmp_gc = data.xc_gain(iC,:);
        if stab(iC)>=0.5
            [pks,loc]=findpeaks(tmp(zero_bin:end),'NPeaks',1);
            [~,mi_loc] = min(tmp); 
            if ~isempty(pks)
            peaks(iC)=loc+zero_bin-1;
            end
            trough(iC)=mi_loc;
            
            [pks,loc]=findpeaks(tmp_gc(zero_bin:end),'NPeaks',1);
            %[~,mi_loc] = min(tmp); 
            if ~isempty(pks)
            peaks_gc(iC)=loc+zero_bin-1;
            end
            
        end
    end
    [~,trough_gc]=min(data.xc_gain,[],2);
    trough_gc(stab<0.5)=nan;
    
    PEAKS_BL = cat(1,PEAKS_BL,nanmedian(peaks));
    PEAKS_GC = cat(1,PEAKS_GC,nanmedian(peaks_gc));
    TROUGHS_BL = cat(1,TROUGHS_BL,nanmedian(trough)); 
    TROUGHS_GC = cat(1,TROUGHS_GC,nanmedian(trough_gc));
    
    XC_BL = cat(1,XC_BL,xc_bl);
    XC_G = cat(1,XC_G,xc_g);
    Cluster_Group = cat(1,Cluster_Group,cg);
    %XC_G2=cat(1,XC_G2,xc_g2);
end
%%
figure
hold on
boundedline(-200:2:200,nanmean(XC_BL),nanstd(XC_BL)/sqrt(size(XC_BL,1)))
boundedline(-200:2:200,nanmean(XC_G),nanstd(XC_G)/sqrt(size(XC_G,1)),'r')
%boundedline(-200:2:200,nanmean(XC_G2),nanstd(XC_G2)/sqrt(size(XC_G2,1)),'g')

legend({'Baseline','0.8 Gain virt cm'})
xlabel('Lag [cm]')
%%
figure
subplot(1,2,1)
idx_stable= Cluster_Group == pca_data.cluster2analyze;
idx_remap = Cluster_Group ~= pca_data.cluster2analyze & ~isnan(Cluster_Group);
dd=TROUGHS_BL-TROUGHS_GC;
histogram(dd(idx_stable),'normalization','probability')
hold on
histogram(dd(idx_remap),'normalization','probability')
title('difference in peaks')
legend('stable','remap')

subplot(1,2,2)
dd=PEAKS_GC-PEAKS_BL;
histogram(dd(idx_stable),'normalization','probability')
hold on
histogram(dd(idx_remap),'normalization','probability')
title('difference in mins')

p_stable = signrank(TROUGHS_BL(idx_stable),TROUGHS_GC(idx_stable));
p_remap = signrank(TROUGHS_BL(idx_remap),TROUGHS_GC(idx_remap));

fprintf('troughs stable %.3f, troughs remap %.3f \n',p_stable,p_remap);

p_stable = signrank(PEAKS_BL(idx_stable),PEAKS_GC(idx_stable));
p_remap = signrank(PEAKS_BL(idx_remap),PEAKS_GC(idx_remap));

fprintf('PEAKS stable %.3f, PEAKS remap %.3f \n',p_stable,p_remap);
fprintf('difference of difference %.3f \n',ranksum(dd(idx_stable),dd(idx_remap)))
%%
figure
hold on
subplot(1,2,1)
idx=Cluster_Group == pca_data.cluster2analyze;
boundedline(-200:2:200,nanmean(XC_BL(idx,:)),nanstd(XC_BL(idx,:))/sqrt(nnz(idx)))

idx=Cluster_Group ~= pca_data.cluster2analyze & ~isnan(Cluster_Group);
boundedline(-200:2:200,nanmean(XC_BL(idx,:)),nanstd(XC_BL(idx,:))/sqrt(nnz(idx)),'r')
%boundedline(-200:2:200,nanmean(XC_G),nanstd(XC_G)/sqrt(size(XC_G,1)),'r')
%boundedline(-200:2:200,nanmean(XC_G2),nanstd(XC_G2)/sqrt(size(XC_G2,1)),'g')

legend({'Baseline Stable','Baseline Remap'})
xlabel('Lag [cm]')

subplot(1,2,2)
idx=Cluster_Group == pca_data.cluster2analyze;
boundedline(-200:2:200,nanmean(XC_G(idx,:)),nanstd(XC_G(idx,:))/sqrt(nnz(idx)))

idx=Cluster_Group ~= pca_data.cluster2analyze & ~isnan(Cluster_Group);
boundedline(-200:2:200,nanmean(XC_G(idx,:)),nanstd(XC_G(idx,:))/sqrt(nnz(idx)),'r')
%boundedline(-200:2:200,nanmean(XC_G),nanstd(XC_G)/sqrt(size(XC_G,1)),'r')
%boundedline(-200:2:200,nanmean(XC_G2),nanstd(XC_G2)/sqrt(size(XC_G2,1)),'g')

legend({'Gain Stable','Baseline Remap'})
xlabel('Lag [cm]')

%plot(-200:2:200,nanmean(XC_G));
%plot(-200:2:200,nanmean(XC_G2));
%%
figure
hold on
xl=[-80 80]
subplot(1,2,1)
idx=Cluster_Group == pca_data.cluster2analyze;
boundedline(-200:2:200,nanmean(XC_BL(idx,:)),nanstd(XC_BL(idx,:))/sqrt(nnz(idx)))
boundedline(-200:2:200,nanmean(XC_G(idx,:)),nanstd(XC_G(idx,:))/sqrt(nnz(idx)),'r')


legend({'Baseline Stable','Gain Stable'})
xlabel('Lag [cm]')
xlabel(sprintf('%.3f',signrank(XC_BL(idx,116),XC_G(idx,116))))

xlim(xl)

subplot(1,2,2)
idx=Cluster_Group ~= pca_data.cluster2analyze & ~isnan(Cluster_Group);
boundedline(-200:2:200,nanmean(XC_BL(idx,:)),nanstd(XC_BL(idx,:))/sqrt(nnz(idx)))
boundedline(-200:2:200,nanmean(XC_G(idx,:)),nanstd(XC_G(idx,:))/sqrt(nnz(idx)),'r')

xlabel(sprintf('%.3f',signrank(XC_BL(idx,116),XC_G(idx,116))))



legend({'Baseline remap','Gain remap'})
%xlabel('Lag [cm]')
xlim(xl)

%plot(-200:2:200,nanmean(XC_G));
%plot(-200:2:200,nanmean(XC_G2));

%%
[~,IDX_BL] =min((XC_BL(:,zero_bin:150)-.2).^2,[],2);
[~,IDX_GC] = min((XC_G(:,zero_bin:150)-.2).^2,[],2);
figure
histogram(IDX_BL)
hold on
histogram(IDX_GC)

figure
idx_stable= Cluster_Group == pca_data.cluster2analyze;
idx_remap = Cluster_Group ~= pca_data.cluster2analyze & ~isnan(Cluster_Group);
subplot(1,2,1)
histogram(IDX_BL(idx_stable))
hold on
histogram(IDX_GC(idx_stable))
legend({'bl','gc'})
title('Width stable')
subplot(1,2,2)
histogram(IDX_BL(idx_remap))
hold on
histogram(IDX_GC(idx_remap))

legend({'bl','gc'})
title('Width remap')

p_stable = signrank(IDX_BL(idx_stable),IDX_GC(idx_stable));
p_remap = signrank(IDX_BL(idx_remap),IDX_GC(idx_remap));

fprintf('width stable %.3f, width remap %.3f \n',p_stable,p_remap);
dd=IDX_BL-IDX_GC;
p = ranksum(dd(idx_stable),dd(idx_remap));

fprintf('difference of difference in width %.3f',p);

figure

plotSpread({dd(idx_stable),dd(idx_remap)})
figure
hold on
histogram(dd(idx_remap),'FaceColor','r','normalization','probability')
histogram(dd(idx_stable),'normalization','probability','FaceColor','b')
legend({'remap','stable'})
set(gca,'XTick',[-4:1:7],'XTickLabel',[-4:1:7]*2)
xlabel('Scale difference [cm]')