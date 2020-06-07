savepath = '/Volumes/Samsung_T5/tbtxcorr_baseline';
pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');

matfiles = dir(fullfile(savepath,'*.mat'));
ops.max_lag = 20;

regions ={'MEC','VISp','RS'};
CorrMat = struct();
SHIFT = struct();
FID = struct();
STAB_BL_Gain = struct();
DEPTH = struct();
SHIFTMat = struct();
FACTORS = struct();
Region = struct();
SessionName = struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    DEPTH.(regions{iR})=[];
    SHIFTMat.(regions{iR})=[];
    FACTORS.(regions{iR}) = [];
    Region.(regions{iR}) = {};
    SessionName.(regions{iR}) = {};
end


idx_m = triu(true(6),1);
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    [~,sn]=fileparts(matfiles(iF).name);
    
    tmp = data_out.corrMat;
    tmp_bl = nan(size(tmp,1),1);
    for iC=1:numel(tmp_bl)
        tmp_m = squeeze(tmp(iC,1:6,1:6));
        tmp_m = tmp_m+diag(nan(6,1));
        tmp_bl(iC)=nanmean(tmp_m(:));
        
    end
    tmp_bl_gain = nanmean(nanmean(tmp(:,1:6,7:10),2),3);
    stab = [tmp_bl';tmp_bl_gain'];
    sm1=data_out.shiftMat; %un altered
    sm2 = data_out.shiftMatS2; %shifted posx, then ran through MGC algorithm
    if isnan(factor(iC))
        sm2(iC,:,:)=nan;
    end

    sm1(isnan(tmp))=nan;
    sm2(isnan(data_out.corrMatS2)) = nan;
    shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
    shift2=nanmean(nanmean(sm2(:,1:6,7:10),2),3);
    shifts = [shift1';shift2'];
    shifts(:,isnan(data_out.factors))=nan;
    tmp_reg = data_out.region;
    XC=data_out.corrMat;
    pm_idx = startsWith(tmp_reg,'VISpm');
    agl_idx = startsWith(tmp_reg,'RSPagl');
    tmp_reg(agl_idx)={'aglRSP'};
    tmp_reg(pm_idx)={'VISm'};
    for iR=1:numel(regions)
        reg_idx = startsWith(tmp_reg,regions{iR});
        if ismember(regions{iR},{'MEC','ECT'})
            mult = -1;
        else
            mult = 1;
        end
        if any(reg_idx)
            CorrMat.(regions{iR}) = cat(1,CorrMat.(regions{iR}),XC(reg_idx,:,:));
            SHIFT.(regions{iR}) = cat(2,SHIFT.(regions{iR}),shifts(:,reg_idx));
            SHIFTMat.(regions{iR}) = cat(1,SHIFTMat.(regions{iR}),sm1(reg_idx,:,:));
            FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
            STAB_BL_Gain.(regions{iR}) = cat(2,STAB_BL_Gain.(regions{iR}),stab(:,reg_idx));
            FACTORS.(regions{iR}) = cat(2,FACTORS.(regions{iR}),data_out.factors(reg_idx));
            Region.(regions{iR}) = cat(2,Region.(regions{iR}),tmp_reg(reg_idx));
            SessionName.(regions{iR}) = cat(2,SessionName.(regions{iR}),repmat({sn},1,nnz(reg_idx)));
        end
    end
    
    
end

%% get cluster1 sessions
ALLC = [];
ALLSHIFTS = [];
ALLR = [];
ALLS = [];
ALLSessionName = {};
for iR=1:3
    sites = FID.(regions{iR});
    usites=unique(sites);
    for iS=1:numel(usites)
        sid = sites==usites(iS);
        stab_idx = STAB_BL_Gain.(regions{iR})(1,:)>.5;
        if nnz(sid & stab_idx)>=5
            idx = sid & stab_idx;
            ALLC=cat(1,ALLC,nanmean(CorrMat.(regions{iR})(idx,:,:)));
            ALLSHIFTS = cat(2,ALLSHIFTS,nanmean(SHIFT.(regions{iR})(:,idx),2));
            ALLR = cat(1,ALLR,iR);
            ALLS = cat(1,ALLS,nanmean(SHIFTMat.(regions{iR})(idx,:,:)));
            ALLSessionName = cat(1,ALLSessionName,unique(SessionName.(regions{iR})(idx)));
        end
    end
end

%%
ALLMEC= ALLC(ALLR==1,:,:);
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
%%

ALLMEC= ALLC(ALLR==1,:,:);
ALLMEC=ALLMEC(:,idx_u);
ALLMEC_ms=ALLMEC-mean(ALLMEC);
clu = zeros(1,size(ALLMEC,1));
score = clu;
for iS=1:size(ALLMEC,1)
    %imagesc(squeeze(ALLMEC(iS,:,:)),[0 .7])
    s = nanmean(nanmean(ALLMEC(iS,:)));
    p = ALLMEC_ms(iS,:);
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


    