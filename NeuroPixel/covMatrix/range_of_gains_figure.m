% savepath = '/Volumes/Samsung_T5/tbtxcorr_contigFR_30cm_data_corrected2';
% 
% matfiles = dir(fullfile(savepath,'*.mat'));
% ops.max_lag = 20;
% 
% regions ={'MEC','VISp','RS'};
% CorrMat = struct();
% SHIFT = struct();
% FID = struct();
% STAB_BL_Gain = struct();
% DEPTH = struct();
% SHIFTMat = struct();
% FACTORS = struct();
% Region = struct();
% SessionName = struct();
% for iR=1:numel(regions)
%     CorrMat.(regions{iR})=[];
%     SHIFT.(regions{iR})=[];
%     STAB_BL_Gain.(regions{iR}) = [];
%     FID.(regions{iR})= [];
%     DEPTH.(regions{iR})=[];
%     SHIFTMat.(regions{iR})=[];
%     FACTORS.(regions{iR}) = [];
%     Region.(regions{iR}) = {};
%     SessionName.(regions{iR}) = {};
% end
% 
% 
% idx_m = triu(true(6),1);
% for iF=1:numel(matfiles)
%     data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
%     [~,sn]=fileparts(matfiles(iF).name);
%     
%     tmp = data_out.corrMat;
%     tmp_bl = nan(size(tmp,1),1);
%     for iC=1:numel(tmp_bl)
%         tmp_m = squeeze(tmp(iC,1:6,1:6));
%         tmp_m = tmp_m+diag(nan(6,1));
%         tmp_bl(iC)=nanmean(tmp_m(:));
%         
%     end
%     tmp_bl_gain = nanmean(nanmean(tmp(:,1:6,7:10),2),3);
%     stab = [tmp_bl';tmp_bl_gain'];
%     sm1=data_out.shiftMat; %un altered
%     sm2 = data_out.shiftMatS2; %shifted posx, then ran through MGC algorithm
%     if isnan(factor(iC))
%         sm2(iC,:,:)=nan;
%     end
%     
%     sm1(isnan(tmp))=nan;
%     sm2(isnan(data_out.corrMatS2)) = nan;
%     shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
%     shift2=nanmean(nanmean(sm2(:,1:6,7:10),2),3);
%     shifts = [shift1';shift2'];
%     shifts(:,isnan(data_out.factors))=nan;
%     tmp_reg = data_out.region;
%     XC=data_out.corrMat;
%     pm_idx = startsWith(tmp_reg,'VISpm');
%     agl_idx = startsWith(tmp_reg,'RSPagl');
%     tmp_reg(agl_idx)={'aglRSP'};
%     tmp_reg(pm_idx)={'VISm'};
%     for iR=1:numel(regions)
%         reg_idx = startsWith(tmp_reg,regions{iR});
%         if ismember(regions{iR},{'MEC','ECT'})
%             mult = -1;
%         else
%             mult = 1;
%         end
%         if any(reg_idx)
%             CorrMat.(regions{iR}) = cat(1,CorrMat.(regions{iR}),XC(reg_idx,:,:));
%             SHIFT.(regions{iR}) = cat(2,SHIFT.(regions{iR}),shifts(:,reg_idx));
%             SHIFTMat.(regions{iR}) = cat(1,SHIFTMat.(regions{iR}),sm1(reg_idx,:,:));
%             FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
%             STAB_BL_Gain.(regions{iR}) = cat(2,STAB_BL_Gain.(regions{iR}),stab(:,reg_idx));
%             FACTORS.(regions{iR}) = cat(2,FACTORS.(regions{iR}),data_out.factors(reg_idx));
%             Region.(regions{iR}) = cat(2,Region.(regions{iR}),tmp_reg(reg_idx));
%             SessionName.(regions{iR}) = cat(2,SessionName.(regions{iR}),repmat({sn},1,nnz(reg_idx)));
%         end
%     end
%     
%     
% end
% 
% %% get cluster1 sessions
% ALLC = [];
% ALLSHIFTS = [];
% ALLR = [];
% ALLS = [];
% ALLSessionName = {};
% for iR=1:3
%     sites = FID.(regions{iR});
%     usites=unique(sites);
%     for iS=1:numel(usites)
%         sid = sites==usites(iS);
%         stab_idx = STAB_BL_Gain.(regions{iR})(1,:)>.5;
%         if nnz(sid & stab_idx)>=5
%             idx = sid & stab_idx;
%             ALLC=cat(1,ALLC,nanmean(CorrMat.(regions{iR})(idx,:,:)));
%             ALLSHIFTS = cat(2,ALLSHIFTS,nanmean(SHIFT.(regions{iR})(:,idx),2));
%             ALLR = cat(1,ALLR,iR);
%             ALLS = cat(1,ALLS,nanmean(SHIFTMat.(regions{iR})(idx,:,:)));
%             ALLSessionName = cat(1,ALLSessionName,unique(SessionName.(regions{iR})(idx)));
%         end
%     end
% end
% 
% tf = triu(true(16),1);
% ALLC(isnan(ALLC))=0;
% [coeff,score,~,~,expl] = pca(ALLC(:,tf));
% 
% % K means clustering
% rng(2); % for reproducibility
% rep_cluster = kmeans(score(:,1:3),3);
% cluster2analyze=3;
% figure
% scatter(score(:,1),score(:,2),15,rep_cluster)
% hold on
% scatter(score(rep_cluster==cluster2analyze,1),score(rep_cluster == cluster2analyze,2),'rx')
% MEC_NAMES = ALLSessionName(ALLR==1);
% MEC_CLUSTERS = rep_cluster(ALLR==1);
% clearvars -except MEC*
%%
matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_range/*.mat')

TG = zeros(212,numel(matfiles));
gains = [0.5, 0.6 0.7 0.8];
SIM_PRE = [];
STAB_PRE = [];
SIM_05 = [];
SID = [];
CORR_MAT=[];
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
    sim_05 = sim_bl_pre;
    onsets_05=strfind(data_out.gain'==0.5,[0 1]) +1;
    corrMats=zeros(nC,10,10,4,2);
    for iG=1:numel(gains)
        
        onsets = strfind(data_out.gain'==gains(iG),[0 1]);
        for iO =1:2
            onset=onsets(iO)+1;
            if gains(iG)==.5
                gain_05_idx = onsets_05(3-iO)+(0:3); %for .5 gain, take other .5 gain
            else
                gain_05_idx = onsets_05(iO)+(0:3); % for other gains, take .5 of same round
            end
            
            bl_pre_idx= onset-6:onset-1;
            sim_bl_pre(:,iG,iO) = nanmean(nanmean(corrMat(:,bl_pre_idx,onset:onset+3),2),3);
            stab_bl_pre(:,iG,iO) =nanmean(nanmean(corrMat(:,bl_pre_idx,bl_pre_idx),2),3);
            sim_05(:,iG,iO) = nanmean(nanmean(corrMat(:,gain_05_idx,onset:onset+3),2),3);
            corrMats(:,:,:,iG,iO)=corrMat(:,(onset-6):(onset+3),(onset-6):(onset+3));
        end
    end
    SIM_PRE = cat(1,SIM_PRE,sim_bl_pre(reg_idx,:,:));
    STAB_PRE = cat(1,STAB_PRE,stab_bl_pre(reg_idx,:,:));
    SIM_05 = cat(1,SIM_05,sim_05(reg_idx,:,:));
    CORR_MAT = cat(1,CORR_MAT,corrMats(reg_idx,:,:,:,:));
    SID = cat(1,SID,iF*ones(nnz(reg_idx),1));
end
%%
figure('Renderer','Painters')
uS = unique(SID);
for iG=1:4
    ALLC=[];
    for iS=1:numel(uS)
        for iO=1:2
            IDX = SID==uS(iS) & squeeze(STAB_PRE(:,iG,iO)>.5);
            if nnz(IDX)>=5
                tmp = squeeze(nanmean(CORR_MAT(IDX,:,:,iG,iO)));
                ALLC=cat(3,ALLC,tmp);
            end
        end
    end
    subplot(1,4,iG)
    hold on
    imagesc(squeeze(nanmean(ALLC,3)),[0 .75])
    axis image
    axis tight
    title(gains(iG))
end
%saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_sim_maps.pdf'))
%%
figure('Position',[440   584   601   214])
subplot(1,2,1)
hold on
uS=unique(SID);
DD=cell(1,4);
for ii=1:4
    DD{ii}=[];
end
for iS=1:numel(uS)
    for iG=1:4
        
        idx1 = STAB_PRE(:,iG,1)>.5 & SID==uS(iS);
        idx2 = STAB_PRE(:,iG,2)>.5 & SID==uS(iS);

        if nnz(idx1)>=5 && nnz(idx2)>5
            p1=nanmean(SIM_05(idx1,iG,1))*.5+nanmean(SIM_05(idx2,iG,2))*.5;
            p2=nanmean(SIM_PRE(idx1,iG,1))*.5 + nanmean(SIM_PRE(idx2,iG,2))*.5;
            DD{iG}=cat(1,DD{iG},p1-p2);
            plot(p1,p2,'.','Color',get_color(gains(iG),100))
        end
    end
end
grid on
axis equal
plot([0 1],[0 1],'k')
xlabel('similarity to .5 map')
ylabel('similarity to bl map')
xlim([.1 .7])
ylim([.1 .7])


subplot(1,2,2)
hold on
for ii=1:4
    plotSpread(DD{ii},'distributionColors',get_color(gains(ii),100),'xValues',ii)
    text(ii-.5,-.2,sprintf('%.3f',signrank(DD{ii},0,'tail','right')))
end
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',gains)
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_blvsgain.pdf'))

%%
matfiles = dir('/Volumes/Samsung_T5/tbtxcorr_range/*.mat');

TG = zeros(212,numel(matfiles));
gains = [0.5, 0.6 0.7 0.8];
SIM_PRE = [];
STAB_PRE = [];
SIM_other_gc = [];
SIM_bl2bl = [];
SID = [];
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    reg_idx = startsWith(data_out.region,'MEC');
    corrMat = data_out.corrMat;
    nC=size(corrMat,1);
    nT=min(numel(data_out.gain),212);
    if nT<150
        continue
    end
    TG(1:nT,iF)=data_out.gain(1:nT);
    sim_bl_pre = zeros(nC,numel(gains),2);
    stab_bl_pre = sim_bl_pre;
    sim_other_gc = sim_bl_pre;
    sim_bl2bl = sim_bl_pre;
    for iG=1:numel(gains)
        
        onsets = strfind(data_out.gain'==gains(iG),[0 1]);
        if numel(onsets)<2
            continue
        end
        for iO = 1:2
            onset_this=onsets(iO)+1;
            onset_other = onsets(3-iO)+1;
            
            
            bl_pre_idx= onset_this-4:onset_this-1;
            
            sim_bl_pre(:,iG,iO) = nanmean(nanmean(corrMat(:,bl_pre_idx,onset_this:onset_this+3),2),3);
            stab_bl_pre(:,iG,iO) =nanmean(nanmean(corrMat(:,bl_pre_idx,bl_pre_idx),2),3);
            
            sim_other_gc(:,iG,iO) = nanmean(nanmean(corrMat(:,onset_other:onset_other+3,onset_this:onset_this+3),2),3);
            sim_bl2bl(:,iG,iO) = nanmean(nanmean(corrMat(:,onset_other-4:onset_other-1,onset_this-4:onset_this-1),2),3);
        end
    end
    SIM_PRE = cat(1,SIM_PRE,sim_bl_pre(reg_idx,:,:));
    STAB_PRE = cat(1,STAB_PRE,stab_bl_pre(reg_idx,:,:));
    SIM_other_gc = cat(1,SIM_other_gc,sim_other_gc(reg_idx,:,:));
    SIM_bl2bl = cat(1,SIM_bl2bl,sim_bl2bl(reg_idx,:,:));
    
    SID = cat(1,SID,iF*ones(nnz(reg_idx),1));
    
end

%%
figure
for ii=1:4
    IDX = STAB_PRE(:,ii,1)>.5;
    p1=squeeze(SIM_other_gc(IDX,ii,1));
    p2=squeeze(SIM_bl2bl(IDX,ii,1));
    subplot(2,4,ii)
    scatter(p1,p2)
    xlabel('Similarity to other gc rep')
    ylabel('Similarity of preceedinb baselines')
    subplot(2,4,ii+4)
    violinplot(p1-p2)
    xlabel(sprintf('%.2f',signrank(p1,p2)))
    ylim([-1 1])
end
%%

%%
figure('Position',[440   584   601   214])

subplot(1,2,1)

hold on
DD=cell(1,4);
for ii=1:4
    DD{ii}=[];
end
for iS=1:numel(uS)
    for ii=1:4
        
        IDX = STAB_PRE(:,ii,1)>.0 & STAB_PRE(:,ii,2)>.0 & SID == uS(iS);
        if nnz(IDX)>=5
            p1=squeeze(SIM_other_gc(IDX,ii,1));
            p2=squeeze(SIM_bl2bl(IDX,ii,1));
            DD{ii}=cat(1,DD{ii},nanmean(p1)-nanmean(p2));
            plot(nanmean(p1),nanmean(p2),'.','Color',get_color(gains(ii),100))
        end
    end
end
xlabel('Similarity to other gc rep')
ylabel('Similarity of preceeding baselines')
grid on
axis equal
xlim([.1 .7])
ylim([.1 .7])
plot([0 1],[0 1],'k')

subplot(1,2,2)
hold on
for ii=1:4
    plotSpread(DD{ii},'distributionColors',get_color(gains(ii),100),'xValues',ii)
    text(ii-.5,-.2,sprintf('%.3f',signrank(DD{ii},0,'tail','right')))
end
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',gains)
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_baselines_gains_noStab_T.pdf'))
