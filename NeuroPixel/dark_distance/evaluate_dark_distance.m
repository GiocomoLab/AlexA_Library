outpath = fullfile(OAK,'distance_tuning_fft');
outpath = '/Users/attialex/distance_tuning_fft';
savepath ='/Volumes/Samsung_T5/attialex/images/dark_distance';
matfiles = dir(fullfile(outpath,'np*.mat'));
hmfig=figure()
DEPTH = [];
TUNED = [];
MOD_DEPTH = [];
PEAK_LOC_XCORR = [];
PEAK_LOC_PXX = [];
for iF = 4%:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    mec_units = startsWith(data_out.region,'MEC');
    depth = data_out.depth(mec_units);
    
    ACG=data_out.ACG(mec_units,:);
    ACG_shuffled = data_out.ACG_shuffled(mec_units,:,:);
    PXX = data_out.PXX(:,mec_units);
    PXX_shuffled = data_out.PXX_shuffled(:,mec_units,:);
    [~,sid] = sort(depth);
    nUnits = nnz(mec_units);
    is_tuned = false(1,nUnits);
    mod_depth = zeros(size(is_tuned));
    peak_loc = zeros(size(is_tuned));
    upper_limits = prctile(ACG_shuffled,99,3);
    upper_limits_pxx = prctile(PXX_shuffled,99,3);
    peak_loc_pxx = peak_loc;
    for iN = 1:nUnits
        [map,mip]=max(PXX(:,iN));
        if upper_limits_pxx(mip,iN)<map
            peak_loc_pxx(iN)=mip;
        end
        
        a=strfind(diff(ACG(iN,:))>0,[0 1]);
        if numel(a)<1
            continue
        end
            
        
        sub = ACG(iN,a(1):end);
        [ma,mi]=max(sub);
        mi = mi+a(1)-1;
        is_tuned(iN) = upper_limits(iN,mi)<ma;
        mod_depth(iN)=ma-ACG(iN,a(1));
        peak_loc(iN)=mi;
%         plot(ACG(iN,:));
%         hold on
%         plot(mi,ACG(iN,mi),'ro')
%         plot(a(1),ACG(iN,a(1)),'rx')
%         if is_tuned(iN)
%             col = 'g';
%         else
%             col = 'r';
%         end
%         plot(upper_limits(iN,:),'Color',col)
% 
%         pause
%         cla
       
    end
    figure(hmfig)
    subplot(1,4,1)
    imagesc(ACG(sid,:))
    
    subplot(1,4,2)
    ACG_tuned = ACG(is_tuned ,:);
    d_tuned = depth(is_tuned );
    [~,sid_t] = sort(d_tuned);
    imagesc(ACG_tuned(sid_t,:))
    
    subplot(1,4,3)
    ACG_tuned = ACG(is_tuned & mod_depth>.2,:);
    d_tuned = depth(is_tuned & mod_depth>.2);
    [~,sid_t] = sort(d_tuned);
    imagesc(ACG_tuned(sid_t,:))
    
    subplot(1,4,4)
    ACG_tuned = ACG(peak_loc_pxx>0,:);
    d_tuned = depth(peak_loc_pxx>0);
    [~,sid_t] = sort(d_tuned);
    imagesc(ACG_tuned(sid_t,:))
    
    saveas(gcf,fullfile(savepath,[sn '.png']))
    clf
    DEPTH = cat(1,DEPTH,depth);
    MOD_DEPTH = cat(1,MOD_DEPTH,mod_depth');
    TUNED = cat(1,TUNED,is_tuned');
    PEAK_LOC_XCORR=cat(1,PEAK_LOC_XCORR,peak_loc');
    PEAK_LOC_PXX = cat(1,PEAK_LOC_PXX,peak_loc_pxx');
end

%%
xbin = 0:ops.SpatialBin:ops.max_lag_autocorr;
bins = 0:6:600;
figure
subplot(3,1,1)
PEAK_LOC_XCORR(isnan(PEAK_LOC_XCORR))=1;
histogram(xbin(PEAK_LOC_XCORR),bins)
xlabel('peak xcorr location')
title('all units')
subplot(3,1,2)
histogram(xbin(PEAK_LOC_XCORR(TUNED==1)),bins)
xlabel('peak xcorr location')
title('tuned units')
subplot(3,1,3)
histogram(xbin(PEAK_LOC_XCORR( MOD_DEPTH>.2)),bins)
xlabel('peak xcorr location')
title('modulated')