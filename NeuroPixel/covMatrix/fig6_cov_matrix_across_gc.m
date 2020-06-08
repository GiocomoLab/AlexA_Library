matfiles = dir('/Volumes/Samsung_T5/attialex/tbtxcorr_range/*.mat');

TG = zeros(212,numel(matfiles));
gains = [0.5, 0.6, 0.7 0.8];
SIM_PRE = [];
STAB_PRE = [];
STAB2_PRE =[];
SIM_05 = [];
SID = [];
CORR_MAT=[];
ACROSS_MAT =[];
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    reg_idx = startsWith(data_out.region,'MEC');
    corrMat = data_out.corrMat;
    nC=size(corrMat,1);
    nT=min(numel(data_out.gain),250);
    if nT<150
        continue
    end
    TG(1:nT,iF)=data_out.gain(1:nT);
    sim_bl_pre = zeros(nC,numel(gains),2);
    stab_bl_pre = sim_bl_pre;
    stab2=stab_bl_pre;
    sim_05 = sim_bl_pre;
    onsets_05=strfind(data_out.gain'==0.5,[0 1]) +1;
    acrossMats = zeros(nC,14,14);
    for iT=1:14
     t1 = onsets_05(1)-11+iT;

        for jT=1:14
        t2=onsets_05(2)-11+jT;
        acrossMats(:,iT,jT)=corrMat(:,t1,t2);
        end
    end
    
    
    corrMats=zeros(nC,20,20,4,2);
    for iG=1:numel(gains)
        
        onsets = strfind(data_out.gain'==gains(iG),[0 1]);
        for iO =1:2
            onset=onsets(iO)+1;
            if gains(iG)==.5
                gain_05_idx = onsets_05(3-iO)+(0:3); %for .5 gain, take other .5 gain
            else
                gain_05_idx = onsets_05(iO)+(0:3); % for other gains, take .5 of same round
            end
            
            bl_pre_idx= onset-4:onset-1;
            sim_bl_pre(:,iG,iO) = nanmean(nanmean(corrMat(:,bl_pre_idx,onset:onset+3),2),3);
            stab_bl_pre(:,iG,iO) =nanmean(nanmean(corrMat(:,bl_pre_idx,bl_pre_idx),2),3);
            stab2(:,iG,iO)=nanmean(nanmean(corrMat(:,onset-6:onset-1,onset-6:onset-1),2),3);
            sim_05(:,iG,iO) = nanmean(nanmean(corrMat(:,gain_05_idx,onset:onset+3),2),3);
            corrMats(:,:,:,iG,iO)=corrMat(:,(onset-10):(onset+9),(onset-10):(onset+9));
        end
    end
    SIM_PRE = cat(1,SIM_PRE,sim_bl_pre(reg_idx,:,:));
    STAB_PRE = cat(1,STAB_PRE,stab_bl_pre(reg_idx,:,:));
    STAB2_PRE = cat(1,STAB2_PRE,stab2(reg_idx,:,:));
    SIM_05 = cat(1,SIM_05,sim_05(reg_idx,:,:));
    CORR_MAT = cat(1,CORR_MAT,corrMats(reg_idx,:,:,:,:));
    ACROSS_MAT = cat(1,ACROSS_MAT,acrossMats(reg_idx,:,:));
    SID = cat(1,SID,iF*ones(nnz(reg_idx),1));
end

%%
%% pre baseline stability - gc stability, only stable cells in baseline
for iG=1
stab_all_bl=[];
stab_all_gc = [];

t_bl = squeeze(nanmean(nanmean(CORR_MAT(:,1:4,1:4,iG,:),2),3));
t_gc = squeeze(nanmean(nanmean(CORR_MAT(:,11:14,11:14,iG,:),2),3));
uID = unique(SID);

for iS=1:numel(uID)
    for iO=1:2
    idx = SID==uID(iS) & squeeze(STAB2_PRE(:,iG,iO)>=.5);
    if nnz(idx)>=5
    t1=nanmean(t_bl(idx,iO));
    t2=nanmean(t_gc(idx,iO));
    stab_all_bl = cat(1,stab_all_bl,t1);
    stab_all_gc = cat(1,stab_all_gc,t2);
    end
    end
end
AV={}
figure('Position',[627   947   613   150])

    %scatter(squeeze(stab_all_bl),squeeze(stab_all_gc))
    AV{2}=squeeze(stab_all_bl);
    AV{1}=squeeze(stab_all_gc);
%     axis image
%     xlim([.2 .65])
%     ylim([.2 .65])
subplot(1,2,1)
plotSpread(AV{1},'xyOri','flipped')
subplot(1,2,2)
plotSpread(AV)
end
set(gca,'XTick',[1,2],'XTickLabel',{'stab gc','diff to pre bl'})
%saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_map_stability_baseline_gain.pdf'))

%%
figure('Position',[680   592   303   506])
plot([1 2],[AV{2}, AV{1}],'.-','MarkerEdgeColor',[0 0 0],'Color',[.5 .5 .5])
xlim([.8 2.2])
set(gca,'XTick',[1 2],'XTickLabel',{'BL','GC'})
xlabel(sprintf('%e',signrank(AV{1},AV{2})))
box off
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_map_stability_baseline_gain.pdf'))


%% pre baseline stability - gc stability, only stable cells in baseline
stab_all_bl=[];
stab_all_gc = [];
t_bl = squeeze(nanmean(nanmean(ACROSS_MAT(:,1:4,1:4 ),2),3));
t_gc = squeeze(nanmean(nanmean(ACROSS_MAT(:,11:14,11:14),2),3));
uID = unique(SID);
for iS=1:numel(uID)
    idx = SID==uID(iS) & squeeze(STAB2_PRE(:,iG,1)>=.5) & squeeze(STAB2_PRE(:,iG,2)>=.5);
    if nnz(idx)>=5
    t1=nanmean(t_bl(idx));
    t2=nanmean(t_gc(idx));
    stab_all_bl = cat(1,stab_all_bl,t1);
    stab_all_gc = cat(1,stab_all_gc,t2);
    end
end
AV={};
figure('Position',[680   592   303   506])
% 
%     subplot(1,2,1)
%     scatter(stab_all_bl,stab_all_gc)
%     
%     xlim([.2 .65])
%     ylim([.2 .65])

plot([1,2],[stab_all_bl,stab_all_gc],'k.')
hold on
plot([1,2],[stab_all_bl,stab_all_gc],'-','Color',[.5 .5 .5])
set(gca,'XTick',[1 2],'XTickLabel',{'BL','GC'})
box off
xlim([.8 2.2])
xlabel(sprintf('%e',signrank(stab_all_gc,stab_all_bl)))
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_map_similarity_acrossStable.pdf'))
%plotSpread({stab_all_gc-stab_all_bl,stab_all_bl,stab_all_gc})
%%
%%  baseline stability - gc stability, all cells in baseline
stab_all_bl=[];
stab_all_gc = [];
t_bl = squeeze(nanmean(nanmean(ACROSS_MAT(:,6:10,6:10 ),2),3));
t_gc = squeeze(nanmean(nanmean(ACROSS_MAT(:,11:14,11:14),2),3));
uID = unique(SID);
for iS=1:numel(uID)
    idx = SID==uID(iS);
    if nnz(idx)>=5
    t1=nanmean(t_bl(idx));
    t2=nanmean(t_gc(idx));
    stab_all_bl = cat(1,stab_all_bl,t1);
    stab_all_gc = cat(1,stab_all_gc,t2);
    end
end
AV={};
figure('Position',[680   592   303   506])
% 
%     subplot(1,2,1)
%     scatter(stab_all_bl,stab_all_gc)
%     
%     xlim([.2 .65])
%     ylim([.2 .65])

plot([1,2],[stab_all_bl,stab_all_gc],'k.')
hold on
plot([1,2],[stab_all_bl,stab_all_gc],'-','Color',[.5 .5 .5])
set(gca,'XTick',[1 2],'XTickLabel',{'BL','GC'})
box off
xlim([.8 2.2])
xlabel(sprintf('%e',signrank(stab_all_gc,stab_all_bl)))
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_map_similarity_acrossAll.pdf'))
%plotSpread({stab_all_gc-stab_all_bl,stab_all_bl,stab_all_gc})
%%
for iG=1
stab_site_bl=[];
stab_site_gc = [];
frac_stable = [];
t_bl = squeeze(nanmean(nanmean(CORR_MAT(:,7:10,7:10,iG,:),2),3));
t_gc = squeeze(nanmean(nanmean(CORR_MAT(:,11:14,11:14,iG,:),2),3));
uID = unique(SID);

for iS=1:numel(uID)
    for iO=1:2
    idx = SID==uID(iS);
    if nnz(idx)>=5
    t1=nanmean(t_bl(idx,iO));
    t2=nanmean(t_gc(idx,iO));
    stab_site_bl = cat(1,stab_site_bl,t1);
    stab_site_gc = cat(1,stab_site_gc,t2);
    n_bl = nnz(t_bl(idx,iO)>=0.5);
    n_gc = nnz(t_gc(idx,iO)>=0.5);
    frac_stable = cat(1,frac_stable,[n_bl, n_gc]./nnz(idx));
    end
    end
end
AV={}
figure
subplot(1,2,1)
scatter(stab_site_bl,stab_site_gc)
axis image
lims = [.1 .7];
xlim(lims)
ylim(lims)
xlabel('stability baseline')
ylabel('stability gain change')
text(0.3,0.3,sprintf('%e',signrank(stab_site_bl,stab_site_gc)))
title('map stability for all neurons')

subplot(1,2,2)
scatter(frac_stable(:,1),frac_stable(:,2))
axis image
lims = [0 .7];
xlim(lims)
ylim(lims)
xlabel('baseline')
ylabel('gain change')
title('frac stable cells')
text(0.3,0.3,sprintf('%e',signrank(frac_stable(:,1),frac_stable(:,2))))

end
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_map_similarity_allCells.pdf'))


%% average similarity maps for blocks that have enough stable cells in both
figure('Renderer','Painters')
uS = unique(SID);
for iG=1:4
    ALLC=nan(10,10,numel(uS));
   
    for iS=1:numel(uS)
        
            IDX = SID==uS(iS) & squeeze(STAB2_PRE(:,iG,1)>=.5) & squeeze(STAB2_PRE(:,iG,2)>=.5);
            if nnz(IDX)>=5
                tmp = squeeze(nanmean(nanmean(CORR_MAT(IDX,5:14,5:14,iG,:),5)));
                ALLC(:,:,iS)=tmp;
            end
      
    end
    subplot(1,4,iG)
    hold on
    imagesc(squeeze(nanmean(ALLC(:,:,:),3)),[0 .75])
    axis image
    axis tight
    title(gains(iG))
    xlabel(nnz(~isnan(ALLC(1,2,:))))
end

%% individual similarity maps, not averaged across blocks in same session
figure('Renderer','Painters')
uS = unique(SID);
for iG=1
    ALLC=nan(16,16,numel(uS),2);
    COMBINED = [];
    for iS=1:numel(uS)
        for iO=1:2
            IDX = SID==uS(iS) & squeeze(STAB2_PRE(:,iG,iO)>.5);
            if nnz(IDX)>=5
                tmp = squeeze(nanmean(CORR_MAT(IDX,5:end,5:end,iG,iO)));
                ALLC(:,:,iS,iO)=tmp;
                COMBINED=cat(1,COMBINED,reshape(tmp,[1 16 16]));
            end
        end
    end
%     subplot(3,4,iG)
%     hold on
%     imagesc(squeeze(nanmean(ALLC(:,:,:,1),3)),[0 .7])
%     axis image
%     axis tight
%     title(gains(iG))
%         subplot(3,4,iG+4)
%     hold on
% 
%         imagesc(squeeze(nanmean(ALLC(:,:,:,2),3)),[0 .7])
%     axis image
%     axis tight
%     title(gains(iG))
    subplot(1,4,iG)
    hold on 
    imagesc(squeeze(nanmean(COMBINED,1)),[0 0.7])
xticks([]);
    yticks([]);
    
    % patches indicating gain value
    trgain = ones(1,16);
    trgain(7:10)=gains(iG);
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
    
    
    
end
%saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_sim_maps.pdf'))
%%
pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');

ALLMEC= COMBINED;
figure
clu = zeros(1,size(ALLMEC,1));
score = clu;
for iS=1:size(ALLMEC,1)
    %imagesc(squeeze(ALLMEC(iS,:,:)),[0 .7])
    s = nanmean(nanmean(ALLMEC(iS,idx_u)));
    p = ALLMEC(iS,idx_u)-pca_data.pca_mu;
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
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_clustering_bar.pdf'))
%% sim matrix across for cells that are stable in both
figure('Renderer','Painters')
uS = unique(SID);
nT=size(ACROSS_MAT,2);
    ALLC=nan(nT,nT,numel(uS));
   cntr=0;
   iG=1;
    for iS=1:numel(uS)
            
            IDX = SID==uS(iS) & squeeze(STAB2_PRE(:,iG,1)>.5) & squeeze(STAB2_PRE(:,iG,2)>.5);
            if nnz(IDX)>=5
                cntr=cntr+1;
                tmp = squeeze(nanmean(ACROSS_MAT(IDX,:,:)));
                ALLC(:,:,iS)=tmp;
            end
        end
  hold on
    imagesc(nanmean(ALLC(5:end,5:end,:),3),[0 .7])
    xlabel(['Block 2 ' num2str(cntr)])
    ylabel('Block 1')
    colorbar
    axis image
    axis tight
    saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_sim_maps_across_6Pre.pdf'))
    %%

%% baseline stability - gc stability
stab_all_bl=[];
stab_all_gc = [];
t_bl = squeeze(nanmean(nanmean(nanmean(CORR_MAT(:,5:10,5:10,[1 4],:),2),3),5));
t_gc = squeeze(nanmean(nanmean(nanmean(CORR_MAT(:,11:14,11:14,[1 4],:),2),3),5));
uID = unique(SID);
for iS=1:numel(uID)
    idx = SID==uID(iS);
    t1=nanmean(t_bl(idx,:,:));
    t2=nanmean(t_gc(idx,:,:));
    stab_all_bl = cat(1,stab_all_bl,t1);
    stab_all_gc = cat(1,stab_all_gc,t2);
end
AV={};
figure('Position',[627   947   613   150])
for ii=1:2
    subplot(1,3,ii)
    scatter(squeeze(stab_all_bl(:,ii,1)),squeeze(stab_all_gc(:,ii,1)))
    AV{ii}=squeeze(stab_all_gc(:,ii,1))-squeeze(stab_all_bl(:,ii,1));
    axis image
    xlim([.2 .65])
    ylim([.2 .65])
end
subplot(1,3,3)
plotSpread(AV)
%% pre baseline stability - gc stability, only stable cells in both baseline
for iG=1
stab_all_bl=[];
stab_all_gc = [];
t_bl = squeeze(nanmean(nanmean(nanmean(CORR_MAT(:,1:4,1:4,iG,:),2),3),5));
t_gc = squeeze(nanmean(nanmean(nanmean(CORR_MAT(:,11:14,11:14,iG,:),2),3),5));
uID = unique(SID);

for iS=1:numel(uID)
    idx = SID==uID(iS) & squeeze(STAB2_PRE(:,iG,1)>=.5) & squeeze(STAB2_PRE(:,iG,2)>=.5);
    if nnz(idx)>=5
    t1=nanmean(t_bl(idx,:,:));
    t2=nanmean(t_gc(idx,:,:));
    stab_all_bl = cat(1,stab_all_bl,t1);
    stab_all_gc = cat(1,stab_all_gc,t2);
    end
end

figure('Position',[627   947   613   150])

    subplot(1,2,1)
    scatter(squeeze(stab_all_bl),squeeze(stab_all_gc))
    AV=squeeze(stab_all_gc)-squeeze(stab_all_bl);
    axis image
    xlim([.2 .65])
    ylim([.2 .65])

subplot(1,2,2)
plotSpread(AV)
end