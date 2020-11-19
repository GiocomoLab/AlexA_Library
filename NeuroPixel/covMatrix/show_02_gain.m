matfiles = dir('/Volumes/T7/attialex/tbtxcorr_02_new/*.mat');
savefig = false;
CORR_MAT=[];
SHIFT_MAT = [];
CORR_MAT_ALL = [];
STAB_ALL = [];
SHIFT_ALL = [];
N_C=[];
ANIMAL = {};
REC = {};
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    corrMat = data_out.corrMat;
    shiftMat = data_out.shiftMat;
    shiftMat(isnan(corrMat))=nan;
    corrMatAll = data_out.corrMat;
    nC=size(corrMat,1);
    reg_idx = startsWith(data_out.region,'VISp');
    
    stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
    stab_tmp = nanmean(nanmean(corrMat(reg_idx,1:6,1:6),2),3);
    tmp_bl_gain = nanmean(nanmean(corrMat(reg_idx,1:6,7:10),2),3);
    stab_all = [stab_tmp tmp_bl_gain];
    shift_all = nanmean(nanmean(shiftMat(reg_idx,1:6,7:10),2),3);
    valid_idx = reg_idx' & stab>=0.5;
    
    if nnz(valid_idx)<5
        continue
    end
    
    CORR_MAT = cat(1,CORR_MAT,(nanmean(corrMat(valid_idx,:,:))));
    SHIFT_MAT = cat(1,SHIFT_MAT,(nanmean(shiftMat(valid_idx,:,:))));
    CORR_MAT_ALL = cat(1,CORR_MAT_ALL,nanmean(corrMatAll(valid_idx,:,:)));
    N_C = cat(1,N_C,[nnz(reg_idx) nnz(valid_idx)]);
    animal = matfiles(iF).name(1:4);
    rec = matfiles(iF).name(1:9);
    ANIMAL = cat(1,ANIMAL,{animal});
    REC = cat(1,REC,{rec});
    STAB_ALL = cat(1,STAB_ALL,stab_all);
    SHIFT_ALL = cat(1, SHIFT_ALL,shift_all);
end
%%
[uA,b,c]=unique(ANIMAL);
[uS,b,c]=unique(REC);

fprintf('%d mice, %d sessions\n',numel(uA),numel(uS))
nSessions = zeros(numel(uS),1);
for iS=1:numel(uS)
    nSessions(iS)=nnz(c==iS);
end
fprintf('%.2f reps per session, max %d \n',mean(nSessions),max(nSessions))
fprintf('%.2f +- %d SEM cells per rep \n',mean(N_C(:,2)),round(std(N_C(:,2))/sqrt(size(N_C,1))))
%%
figure
scatter(STAB_ALL(:,1),STAB_ALL(:,2),55,SHIFT_ALL,'.')
xlabel('stab bl pre')
ylabel('stab gain change')
set(gca,'CLim',[-5 5])
axis image
xlim([-.2 1])
ylim([-.2 1])
%% post baseline stability - gc stability, only stable cells in baseline


t_bl = squeeze(nanmean(nanmean(CORR_MAT_ALL(:,1:4,1:4),2),3));
t_gc = squeeze(nanmean(nanmean(CORR_MAT_ALL(:,11:14,11:14),2),3));

AV={}
figure('Position',[627   947   613   150])

%scatter(squeeze(stab_all_bl),squeeze(stab_all_gc))
AV{2}=squeeze(t_bl);
AV{1}=squeeze(t_gc);

figure('Position',[680   592   303   506])      
plot([1 2],[AV{2}, AV{1}],'.-','MarkerEdgeColor',[0 0 0],'Color',[.5 .5 .5])
xlim([.8 2.2])
set(gca,'XTick',[1 2],'XTickLabel',{'BL','GC'})
xlabel(sprintf('%e',signrank(AV{1},AV{2})))
box off
if savefig
    saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_map_stability_baseline_gain.pdf'))
end
%% similarity map
figure
hold on
imagesc(squeeze(nanmean(CORR_MAT,1)),[0 0.7])
xticks([]);
yticks([]);

% patches indicating gain value
trgain = ones(1,16);
trgain(7:10)=0.2;
num_tr = 16;
opt.num_tr_bl=6;
opt.num_tr_gc = 4;
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
        'EdgeColor','none');
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
        'EdgeColor','none');
end
xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);

% white lines delineated gain change trials
for tr = 0.5+[opt.num_tr_bl opt.num_tr_bl+opt.num_tr_gc]
    plot([0.5 num_tr+0.5],[tr tr],'w-');
    plot([tr tr],[0.5 num_tr+0.5],'w-');
end

xlim([-1 num_tr+0.5]);
ylim([-1 num_tr+0.5]);
axis square;
axis off
if savefig
    saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_sim_mapsRS.pdf'))
end
%%
%% shift map
figure
hold on
imagesc(squeeze(nanmean(SHIFT_MAT,1)),[-5 5])
xticks([]);
yticks([]);

% patches indicating gain value
trgain = ones(1,16);
trgain(7:10)=0.5;
num_tr = 16;
opt.num_tr_bl=6;
opt.num_tr_gc = 4;
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
        'EdgeColor','none');
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
        'EdgeColor','none');
end
xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);

% white lines delineated gain change trials
for tr = 0.5+[opt.num_tr_bl opt.num_tr_bl+opt.num_tr_gc]
    plot([0.5 num_tr+0.5],[tr tr],'w-');
    plot([tr tr],[0.5 num_tr+0.5],'w-');
end

xlim([-1 num_tr+0.5]);
ylim([-1 num_tr+0.5]);
axis square;
axis off
colorbar
title('RS')
if savefig
    saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_shift_mapsRS.pdf'))
end
%%
%pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/rep_clusters.mat');
pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');

ALLMEC= CORR_MAT;
figure
clu = zeros(1,size(ALLMEC,1));
score = clu;
idx_u = triu(true(16),1)
for iS=1:size(ALLMEC,1)
    %imagesc(squeeze(ALLMEC(iS,:,:)),[0 .7])
    s = nanmean(nanmean(ALLMEC(iS,idx_u)));
    p = ALLMEC(iS,idx_u)-mu;%mpca_data.pca_mu;
    c = p*pca_data.pca_coeff(:,1:3);
    d = vecnorm(c-pca_data.centroids,2,2);
    [~,clu(iS)]=min(d);
    score(iS)=s;
    %xlabel(sprintf('clu = %d, %.2f',clu(iS),s))
    %pause
    %cla
end
idx_1=clu==3;
idx_3=clu==1;
clu(idx_1)=1;
clu(idx_3)=3;
figure
counts =histcounts(clu,3);
barh([1 2],[counts;counts],'stacked')
xlabel(sprintf('%.2f in clu 1, %.2f abouve .55',counts(1)/sum(counts),nnz(score>.55)/numel(score)))
if savefig
    saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_clustering_bar.pdf'))
end