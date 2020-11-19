savepath = '/Volumes/T7/attialex/tbtxcorr_0.5_withPB';

matfiles = dir(fullfile(savepath,'*.mat'));
ops.max_lag = 20;

regions ={'MEC','VISp','RS'};
CorrMat = struct();
SHIFT = struct();
FID = struct();
STAB_BL_Gain = struct();
DEPTH = struct();
SHIFTMat = struct();
STAB_PB = struct();
Region = struct();
SessionName = struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    DEPTH.(regions{iR})=[];
    SHIFTMat.(regions{iR})=[];
    STAB_PB.(regions{iR}) = [];
    Region.(regions{iR}) = {};
    SessionName.(regions{iR}) = {};
end


idx_m = triu(true(6),1);
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    [~,sn]=fileparts(matfiles(iF).name);
    stab_pb = nanmean(nanmean(data_out.corrMatPB(:,1:20,1:20),2),3);
    stab_bl = nanmean(nanmean(data_out.corrMatBL(:,1:20,1:20),2),3);
    tmp = data_out.corrMat;
    tmp_bl = nan(size(tmp,1),1);
    for iC=1:numel(tmp_bl)
        tmp_m = squeeze(tmp(iC,1:6,1:6));
        tmp_m = tmp_m+diag(nan(6,1));
        tmp_bl(iC)=nanmean(tmp_m(:));
        
    end
    tmp_bl_gain = nanmean(nanmean(tmp(:,1:6,7:10),2),3);
    stab = [tmp_bl';tmp_bl_gain';stab_bl'];
    sm1=data_out.shiftMat; %un altered
    

    sm1(isnan(tmp))=nan;
    shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
    shifts = [shift1'];
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
            STAB_PB.(regions{iR}) = cat(1,STAB_PB.(regions{iR}),stab_pb(reg_idx));
            FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
            STAB_BL_Gain.(regions{iR}) = cat(2,STAB_BL_Gain.(regions{iR}),stab(:,reg_idx));
            Region.(regions{iR}) = cat(2,Region.(regions{iR}),tmp_reg(reg_idx));
            SessionName.(regions{iR}) = cat(2,SessionName.(regions{iR}),repmat({sn},1,nnz(reg_idx)));
        end
    end
    
    
end

%%
reg ={'VISp','RS'}
for iR=1:2
figure
subplot(1,4,1)
scatter(STAB_BL_Gain.(reg{iR})(1,:),STAB_PB.(reg{iR}),35,SHIFT.(reg{iR}),'.')
xlabel('stability BL Pre gain session')
ylabel('stability playback')
axis image
xlim([-.2 1])
ylim([-.2 1])
set(gca,'CLim',[-5 5])
title(reg{iR})


subplot(1,4,2)
scatter(STAB_BL_Gain.(reg{iR})(3,:),STAB_PB.(reg{iR}),35,SHIFT.(reg{iR}),'.')
xlabel('stability BL')
ylabel('stability playback')
axis image
xlim([-.2 1])
ylim([-.2 1])
set(gca,'CLim',[-5 5])

subplot(1,4,3)
scatter(STAB_BL_Gain.(reg{iR})(2,:),STAB_PB.(reg{iR}),35,SHIFT.(reg{iR}),'.')
xlabel('stability Gain')
ylabel('stability playback')
axis image
xlim([-.2 1])
ylim([-.2 1])
set(gca,'CLim',[-5 5])

subplot(1,4,4)
scatter(SHIFT.(reg{iR}),STAB_PB.(reg{iR}),35,STAB_BL_Gain.(reg{iR})(1,:),'.')
xlabel('shift gain')
ylabel('stability PB')
grid on
end