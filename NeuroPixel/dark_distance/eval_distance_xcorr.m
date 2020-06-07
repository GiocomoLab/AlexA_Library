OAK = '/Volumes/Samsung_T5/attialex';
outpath = fullfile(OAK,'distance_tuning_xcorr_only');
%outpath = '/Volumes/Samsung_T5/attialex/';
%savepath ='/Volumes/Samsung_T5/attialex/images/dark_distance';
%savepath =fullfile(OAK,'images_dark_distance_welch');
if ~isfolder(savepath)
    mkdir(savepath)
end
matfiles = dir(fullfile(outpath,'np*.mat'));
hmfig=figure();
DEPTH = [];
TUNED = [];
MOD_DEPTH = [];
PEAK_LOC_XCORR = [];
PEAK_LOC_PXX = [];
SID = [];
FRAC = [];
REGION = {};
for iF = 1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    mec_units = startsWith(data_out.region,{'MEC','ECT'});
    tip_distance = data_out.depth(mec_units);
    if ~isnan(data_out.mec_entry)
    depth = tip_distance - data_out.mec_entry;
    else
        depth = tip_distance - nan; % we dont have a z2 for some cases, could estimate it in the future
    end
    
    ACG=data_out.ACG(mec_units,:);
    
    [~,sid] = sort(depth);
    nUnits = nnz(mec_units);
    is_tuned = false(1,nUnits);
    mod_depth = zeros(size(is_tuned));
    peak_loc = zeros(size(is_tuned));
    upper_limits = squeeze(data_out.acg_upper_limits(:,:,end)); 
    
    for iN = 1:nUnits
        
        
        a=strfind(diff(ACG(iN,:))>0,[0 1]);
        if numel(a)<1 %no elbow point detected, skip
            continue
        end
        
        if a(1)>250 % elbow point beyond 500cm, probably noise, skip
            continue
        end
            
        
        sub = ACG(iN,a(1):250);
        [ma,mi]=max(sub);
        mi = mi+a(1)-1;
        is_tuned(iN) = upper_limits(iN,mi)<ma;
        mod_depth(iN)=ma-ACG(iN,a(1));
        peak_loc(iN)=mi;

       
    end
%     figure(hmfig)
%     subplot(1,3,1)
%     imagesc(ACG(sid,:))
%     
%     subplot(1,3,2) 
%     ACG_tuned = ACG(is_tuned ,:);
%     d_tuned = depth(is_tuned );
%     [~,sid_t] = sort(d_tuned);
%     imagesc(ACG_tuned(sid_t,:))
%     hold on
%     %plot(peak_loc(is_tuned),sid_t,'r*')
%     peaks_tuned = peak_loc(is_tuned);
% y_vec = 1:nnz(is_tuned);
% plot(peaks_tuned(sid_t),y_vec,'ro')
% 
%     subplot(1,3,3)
%     ACG_tuned = ACG(is_tuned & mod_depth>.1,:);
%     d_tuned = depth(is_tuned & mod_depth>.1);
%     [~,sid_t] = sort(d_tuned);
%     imagesc(ACG_tuned(sid_t,:))
%     hold on
%         peaks_tuned = peak_loc(is_tuned & mod_depth>.1);
% y_vec = 1:nnz(is_tuned & mod_depth>.1);
% plot(peaks_tuned(sid_t),y_vec,'ro')
% 
%     
%     
%     %saveas(gcf,fullfile(savepath,[sn '.png']))
%     pause
%     clf
    DEPTH = cat(1,DEPTH,depth);
    MOD_DEPTH = cat(1,MOD_DEPTH,mod_depth');
    TUNED = cat(1,TUNED,is_tuned');
    PEAK_LOC_XCORR=cat(1,PEAK_LOC_XCORR,peak_loc');
    SID = cat(1,SID,iF*ones(numel(mod_depth),1));
    REGION = cat(1,REGION,data_out.region');
end



%% 2d scatter plot and histogram of peak locatoions, for all cells
IDX = ( MOD_DEPTH>.1 & TUNED == 1);

figure
scatter(xbin(PEAK_LOC_XCORR(IDX)),DEPTH(IDX))

figure
subplot(1,6,[1:5])
h=histogram2(xbin(PEAK_LOC_XCORR(IDX)),DEPTH(IDX)',[70 70],'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
subplot(1,6,6)
aid = unique(SID);
midpoints =  -2300:200:1500;
half_width = 100;
min_number_cells = 5;
xlabel('Peak Location')
ylabel('Depth relative to MEC entry')
fraction_tuned = nan(numel(aid));
allSitesFraction = nan(numel(aid),numel(midpoints));
for ii=1:numel(aid)
    site_idx = SID==ii;
    valid_idx = IDX & site_idx;
    vals = nan(1,numel(midpoints));
    valsF = nan(1,numel(midpoints));
    for iv=1:numel(vals)
    idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    depth_idx  = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width;
    if sum(idx)>min_number_cells
    valsF(iv)=nnz(idx)/nnz(depth_idx&site_idx);
    end
    end
    %plot(vals,edges(1:end-1))
    allSites(ii,:)=vals;
    allSitesFraction(ii,:)=valsF;
end
errF = nanstd(allSitesFraction)./sqrt(sum(~isnan(allSitesFraction)));
errorbar(nanmean(allSitesFraction),midpoints,[],[],errF,errF,'r')
xlabel('Frac distance tuned')
xlim([0 1])
ylim(h.YBinLimits)
set(gca,'YTickLabel','')
%% histogram of peak locations, for MEC cells only
MEC_idx = startsWith(REGION,'MEC');
xbin = 0:ops.SpatialBin:ops.max_lag_autocorr;
bins = 0:6:ops.max_lag_autocorr;
figure
subplot(3,1,1)
PEAK_LOC_XCORR((PEAK_LOC_XCORR)==0)=1;
histogram(xbin(PEAK_LOC_XCORR(MEC_idx)),bins)
xlabel('peak xcorr location')
title('all units')
subplot(3,1,2)
histogram(xbin(PEAK_LOC_XCORR(TUNED==1 & MEC_idx)),bins)
xlabel('peak xcorr location')
title('tuned units')
subplot(3,1,3)
histogram(xbin(PEAK_LOC_XCORR( MOD_DEPTH>.1 & TUNED == 1 & MEC_idx)),bins)
xlabel('peak xcorr location')
title('modulated')

hold on; for ii=1:7; xline(75*sqrt(2)^(ii-2));end

%% fraction of distance tuned units in MEC
figure
IDX = MOD_DEPTH>.1 & TUNED == 1 & MEC_idx;
aid = unique(SID);
frac = nan(numel(aid),1);
for iS=1:numel(aid)
    frac(iS)=nnz(IDX & SID==aid(iS))/nnz(SID==aid(iS));
end
plotSpread(frac,'SpreadWidth',0.2)

boxplot(frac)
title('Based on XCorr')