matfiles = dir('/Volumes/Samsung_T5/attialex/autocorr2/*.mat');

%%
XC_BL=[];
XC_G=[];
XC_G2=[];
x_vec = -200:2:200;
zero_bin = find(x_vec==0);
PEAKS = [];
for iF=1:numel(matfiles)
    data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    stab=nanmean(nanmean(data.corrMat(:,1:6,1:6),2),3);
    stab_idx=stab>=0.5;
    if nnz(stab_idx)<5
        continue
    end
    xc_bl = nanmean(data.xc_bl(stab_idx,:));
    xc_g = nanmean(data.xc_gain(stab_idx,:));
    xc_g2 = nanmean(data.xc_gain_corrected(stab_idx,:));
    nC=numel(stab);
    peaks = nan(nC,1);
    for iC=1:nC
        tmp = data.xc_bl(iC,:);
        if stab(iC)>=0.5
            [pks,loc]=findpeaks(tmp(zero_bin:end),'NPeaks',1);
            peaks(iC)=loc+zero_bin-1;
            
        end
        
        XC_BL = cat(1,XC_BL,xc_bl);
        XC_G = cat(1,XC_G,xc_g);
        XC_G2=cat(1,XC_G2,xc_g2);
        PEAKS = cat(1,PEAKS,median(peaks));
    end
    %%
    figure
    hold on
    boundedline(-200:2:200,nanmean(XC_BL),nanstd(XC_BL)/sqrt(size(XC_BL,1)))
    boundedline(-200:2:200,nanmean(XC_G),nanstd(XC_G)/sqrt(size(XC_G,1)),'r')
    boundedline(-200:2:200,nanmean(XC_G2),nanstd(XC_G2)/sqrt(size(XC_G2,1)),'g')
    
    legend({'Baseline','0.5 Gain virt cm','0.5 Gain real cm'})
    xlabel('Lag [cm]')
    %plot(-200:2:200,nanmean(XC_G));
    %plot(-200:2:200,nanmean(XC_G2));
    %% tuned vs non tuned across cells
    XC_BLTuned=[];
    XC_BLNonTuned = [];
    XC_GTuned=[];
    XC_GNonTuned=[];
    XC_G2Tuned=[];
    XC_G2NonTuned = [];
    
    for iF=1:numel(matfiles)
        data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
        if any(isnan(data.is_tuned))
            continue
        end
        
        stab=nanmean(nanmean(data.corrMat(:,1:6,1:6),2),3);
        stab_idx=stab>=0.5;
        if nnz(stab_idx)<5
            continue
        end
        tuned_idx = data.is_tuned==1;
        %t1 = nanmean(data.xc_bl(stab_idx & tuned_idx,:));
        %t2 = nanmean(data.xc_bl(stab_idx & ~tuned_idx,:));
        t1 = (data.xc_bl(stab_idx & tuned_idx,:));
        t2 = (data.xc_bl(stab_idx & ~tuned_idx,:));
        XC_BLTuned = cat(1,XC_BLTuned,t1);
        XC_BLNonTuned = cat(1,XC_BLNonTuned,t2);
        
        t1 = (data.xc_gain(stab_idx & tuned_idx,:));
        t2 = (data.xc_gain(stab_idx & ~tuned_idx,:));
        XC_GTuned = cat(1,XC_GTuned,t1);
        XC_GNonTuned = cat(1,XC_GNonTuned,t2);
        
        
    end
    
    
end
%%
figure
subplot(1,2,1)
plot(-200:2:200,nanmean(XC_BLTuned))
hold on
plot(-200:2:200,nanmean(XC_BLNonTuned))
legend({'Tuned','Not Tuned'})
xlabel('lag')
title('Autocorr Baseline')
subplot(1,2,2)
plot(-200:2:200,nanmean(XC_GTuned))
hold on
plot(-200:2:200,nanmean(XC_GNonTuned))
legend({'Tuned','Not Tuned'})
xlabel('lag')
title('Autocorr Gain')
%%
figure
subplot(1,2,1)
plot(-200:2:200,nanmean(XC_BLTuned))
hold on
plot(-200:2:200,nanmean(XC_GTuned))
legend({'Lag BL Tuned','Lag Gain Tuned'})
xlabel('cm')
title('Distance tuned')
subplot(1,2,2)
plot(-200:2:200,nanmean(XC_BLNonTuned))
hold on
plot(-200:2:200,nanmean(XC_GNonTuned))
legend({'Lag BL NonTuned','Lag Gain Non Tuned'})
xlabel('cm')
title('Not Distance Tuned')
