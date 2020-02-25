files = dir('Z:\giocomo\attialex\distance_coding3\data/np*dark*.mat');
idx = true(1,numel(files));
for iF=1:numel(files)
    if contains(files(iF).name,'playback') ||contains(files(iF).name,'reward')
        idx(iF)=false;
    end
end
files = files(idx);
        
PEAKS=[];
PEAKValue=[];
QUANT = [];
DEPTH = [];
PXXPeak = [];
PXXSort = [];
AMP=[];
centers = linspace(1000,-2000,30);
ACG_LIST = {numel(centers,1)};
SID=[];
SMALLMOD = zeros(numel(files));
THETA = [];
THETA_NORM = [];
THETA_NORM2= [];
FIRING_RATE = [];
AVG_CORR = [];
for iL = 1:numel(centers)
    ACG_LIST{iL}=[];
end
for iF=1:numel(files)
    clear firing_rate
    clear theta
    clear avg_corr
    load(fullfile(files(iF).folder,files(iF).name))
    [a,b]=sort(mec_depth);
    %     imagesc(ACG(b,:),[0 0.4])
    %     [x,y]=ginput()
    %     SMALLMOD(iF)=x(end);
    v_idx =1./f_vec > 25 &  1./f_vec<750;
    %PXX(upper_bound_pxx'>PXX)=0;
    [ma,mi]=max(PXX(v_idx,:));
    offset = strfind(v_idx,[0 1]);
    mi = mi+offset;
    
    [pval,psid]=sort(PXX(v_idx,:),'descend');
    psid = psid+offset;
    tmp = zeros(numel(mi),3,3);
    for iC=1:numel(mi)
        tmp(iC,:,1)=pval(1:3,iC);
        tmp(iC,:,2)=psid(1:3,iC);
        tmp(iC,:,3)=upper_bound_pxx(iC,psid(1:3,iC));
    end
    
%     tmp = zeros(numel(mi),3,3);
%     i_vec=-1:1;
%     for iC=1:numel(mi)
%         tmp(iC,:,1)=PXX(mi(iC)+i_vec,iC);
%         tmp(iC,:,2)=[mi(iC) mi(iC) mi(iC)];
%         tmp(iC,:,3)=upper_bound_pxx(iC,i_vec+mi(iC));
%     end
    try
    AVG_CORR = cat(1,AVG_CORR,avg_corr);
    catch
        sprintf('no corr for %d \n',iF)
    end
    PXXSort=cat(1,PXXSort,tmp);
    THETA = cat(2,THETA,theta);
    THETA_NORM = cat(2,THETA_NORM,normalized_theta);
    THETA_NORM2 = cat(2,THETA_NORM2,normalized_theta2);
    FIRING_RATE = cat(2,FIRING_RATE,firing_rate');
    for iC=1:numel(peak_list)
        miACG = min(ACG(iC,:));
        
        SID(end+1)=iF;
        if ~isempty(peak_list(iC).peak_loc)
            [map,mip]=max(peak_list(iC).peak_val);
            %mip=1;
            PEAKS(end+1)=peak_list(iC).peak_loc(mip);
            PEAKValue(end+1)=peak_list(iC).peak_val(mip);
            QUANT(end+1)=peak_list(iC).quantile(mip);
            AMP(end+1)=map-miACG;
            DEPTH(end+1)=mec_depth(iC);
            
            [~,icenter] =min(abs(mec_depth(iC)-centers));
            if peak_list(iC).peak_val(1)>peak_list(iC).quantile(1)
                ACG_LIST{icenter}=cat(1,ACG_LIST{icenter},ACG(iC,:));
            end
        else
            PEAKS(end+1)=nan;
            PEAKValue(end+1)=nan;
            QUANT(end+1)=nan;
            DEPTH(end+1)=nan;
            AMP(end+1)=nan;
            
        end
        
        if ma(iC)>upper_bound_pxx(iC,mi(iC))
            PXXPeak(end+1)=1/f_vec(mi(iC));
        else
            PXXPeak(end+1)=nan;
        end
    end
    
end
% %% normalize for xcorr peaks
% aid = unique(SID);
% size_idx=PEAKValue>(QUANT+.00);
% frac =.2;
% peaks_norm = [];
% figure
% for ii=1:numel(aid)
%     idx = SID==ii;
%     valid_idx = size_idx & idx;
%     depth_this = DEPTH(valid_idx);
%     peaks_this = PEAKS(valid_idx);
%     [~,sid]=sort(depth_this,'descend');
%     n=round(numel(sid)*frac);
%     val = median(peaks_this(sid(1:n)));
%     peaks_norm=cat(2,peaks_norm,peaks_this/val);
%     %     subplot(1,2,1)
%     %     histogram(peaks_this)
%     %     subplot(1,2,2)
%     %     histogram(peaks_this/val,20);
%     %     pause
%     %     clf
% end
% figure
% scatter(peaks_norm,DEPTH(size_idx))
% figure
% histogram(peaks_norm,100)
% figure
% 
% h = histogram2(peaks_norm,DEPTH(size_idx),[12 12],'FaceColor','flat');

% %% normalize for pxx peaks
% aid = unique(SID);
% size_idx=~isnan(PXXPeak);
% frac =.2;
% peaks_norm = [];
% figure
% for ii=1:numel(aid)
%     idx = SID==ii;
%     valid_idx = size_idx & idx;
%     depth_this = DEPTH(valid_idx);
%     peaks_this = PXXPeak(valid_idx);
%     %[~,sid]=sort(depth_this,'descend');
%     %n=round(numel(sid)*frac);
%     %val = median(peaks_this(sid(1:n)));
%     [~,sid]=sort(peaks_this);
%     n=round(frac*numel(sid));
%     val = mean(peaks_this(sid(1:n)));
%     peaks_norm=cat(2,peaks_norm,peaks_this/val);
%     %     subplot(1,2,1)
%     %     histogram(peaks_this)
%     %     subplot(1,2,2)
%     %     histogram(peaks_this/val,30);
%     %     pause
%     %     clf
% end
% figure
% scatter(peaks_norm,DEPTH(size_idx))
% figure
% histogram(peaks_norm,100)
% figure
% 
% h = histogram2(peaks_norm,DEPTH(size_idx),[12 12],'FaceColor','flat');
FR_THRESHOLD=1;
%%
figure('Position',[34         638        1508         340])
subplot(1,3,1)
histogram(PEAKS(PEAKValue>(QUANT+.00) & AMP>0 & FIRING_RATE>FR_THRESHOLD)*5,151)
title('Xcorr Based')
subplot(1,3,2)
h=histogram(PXXPeak(FIRING_RATE>FR_THRESHOLD),111,'BinLimits',[35 580.3000]);
title('PXX Based')
subplot(1,3,3)
size_idx = all(PXXSort(:,:,1)>PXXSort(:,:,3),2) & FIRING_RATE'>FR_THRESHOLD;

h=histogram(PXXPeak(size_idx),111,'BinLimits',[35 800.3000]);
title('PXX Based')


%%
% figure
% valididx = ~isnan(PXXPeak) & ~isnan(DEPTH);
% scatplot(PXXPeak(valididx),DEPTH(valididx))

%%
figure('Position', [249         482        1484         391])
subplot(1,3,1)
size_idx = PEAKValue>(QUANT+.0) & AMP>0.0 & FIRING_RATE>FR_THRESHOLD;
h = histogram2(PEAKS(size_idx)*5,DEPTH(size_idx),[70 70],'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
set(gca,'CLim',[1 15])
xlabel('Peak location [cm]')
ylabel('Distance from MEC entry [um]')
colorbar
set(gcf,'Render','Painters')

axis square
subplot(1,3,2)
imagesc(flipud(h.Values'),[1 15])
colorbar
axis image
subplot(1,3,3)

aid = unique(SID);
frac_dist_coding_xcorr=zeros(1,numel(aid));
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    frac_dist_coding_xcorr(ii)=sum(valid_idx)/sum(idx);
end
plotSpread(frac_dist_coding_xcorr','SpreadWidth',0.2)

boxplot(frac_dist_coding_xcorr)
title('Based on XCorr')

%%

figure('Position', [249         482        1484         391])
subplot(1,3,1)
size_idx = ~isnan(PXXPeak) & FIRING_RATE>FR_THRESHOLD;
h = histogram2(PXXPeak(size_idx),DEPTH(size_idx),[70 70],'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
set(gca,'CLim',[1 15])
xlabel('Peak location [cm]')
ylabel('Distance from MEC entry [um]')
colorbar
set(gcf,'Render','Painters')

axis square
subplot(1,3,2)
imagesc(flipud(h.Values'),[1 15])
colorbar
axis image
subplot(1,3,3)

aid = unique(SID);
frac_dist_coding=zeros(1,numel(aid));
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    frac_dist_coding(ii)=sum(valid_idx)/sum(idx);
end
plotSpread(frac_dist_coding','SpreadWidth',0.2)

boxplot(frac_dist_coding)
title('Based on PSD')
%%

figure('Position', [249         482        1484         391])
subplot(1,3,1)
%size_idx = PXXSort(:,1,1)>PXXSort(:,1,3) & FIRING_RATE'>FR_THRESHOLD;
size_idx = all(PXXSort(:,:,1)>PXXSort(:,:,3),2) & FIRING_RATE'>FR_THRESHOLD;
h = histogram2(1./f_vec(PXXSort(size_idx,1,2)),DEPTH(size_idx),[70 70],'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
set(gca,'CLim',[0 10])
xlabel('Peak location [cm]')
ylabel('Distance from MEC entry [um]')
colorbar
set(gcf,'Render','Painters')

axis square
subplot(1,3,2)
imagesc(flipud(h.Values'),[1 15])
colorbar
axis image
subplot(1,3,3)

aid = unique(SID);
frac_dist_coding2=zeros(1,numel(aid));
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx';
    frac_dist_coding2(ii)=sum(valid_idx)/sum(idx);
end
plotSpread(frac_dist_coding2','SpreadWidth',0.2)

boxplot(frac_dist_coding2)
title('Based on PSD')

%%
aid = unique(SID);
midpoints= -2300:200:1500;
half_width = 100;
idx_value = ~isnan(PXXPeak)& FIRING_RATE>1;
%idx_value = all(PXXSort(:,:,1)>PXXSort(:,:,3),2) & FIRING_RATE'>FR_THRESHOLD;
%idx_value = idx_value';
%idx_value = true(size(PXXPeak)) & FIRING_RATE>FR_THRESHOLD;
%idx_value = PEAKValue>QUANT & AMP>.2;



allSitesRaw = nan(numel(aid),numel(midpoints));
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = idx_value & idx;
    vals = nan(1,numel(midpoints));
    for iv=1:numel(vals)
    idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    if sum(idx)>5
    vals(iv)=nanmean(THETA(idx));
    end
    end
    %plot(vals,edges(1:end-1))
    allSitesRaw(ii,:)=vals;
end




figure
subplot(1,3,1)
%plot(nanmean(allSites),edges(1:end-1))
err = nanstd(allSitesRaw)./sqrt(sum(~isnan(allSitesRaw)));
errorbar(nanmean(allSitesRaw),midpoints,[],[],err,err)
ylabel('Depth')
xlabel('Theta Power')




allSites = nan(numel(aid),numel(midpoints));
allSitesN = allSites;
allSitesFraction = allSites;
for ii=1:numel(aid)
    site_idx = SID==ii;
    valid_idx = idx_value & site_idx;
    vals = nan(1,numel(midpoints));
    valsF = nan(1,numel(midpoints));
    for iv=1:numel(vals)
    idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    depth_idx  = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width;
    if sum(idx)>5
    vals(iv)=nanmean(THETA_NORM(idx));
    valsF(iv)=nnz(idx)/nnz(depth_idx&site_idx);
    end
    allSitesN(ii,iv)=sum(idx);
    end
    %plot(vals,edges(1:end-1))
    allSites(ii,:)=vals;
    allSitesFraction(ii,:)=valsF;
end
%set(gca,'XScale','log')

ax1=subplot(1,3,2);

%err = nanstd(allSites)./sqrt(numel(aid));
enough_idx = sum(~isnan(allSites));
enough_idx = enough_idx>0;
err = nanstd(allSites)./sqrt(sum(~isnan(allSites)));

errorbar(nanmean(allSites(:,enough_idx)),midpoints(enough_idx),[],[],err(enough_idx),err(enough_idx))
ylabel('Depth')
xlabel('Normalized theta Power')

subplot(1,3,3)
errF = nanstd(allSitesFraction)./sqrt(sum(~isnan(allSitesFraction)));
errorbar(nanmean(allSitesFraction(:,enough_idx)),midpoints(:,enough_idx),[],[],errF(:,enough_idx),errF(:,enough_idx),'r')
ylabel('Depth')
xlabel('Frac distance tuned')



%% firing rate and such

aid = unique(SID);
midpoints= -2300:200:1500;
half_width = 100;
allSitesRate = nan(numel(aid),numel(midpoints));

for ii=1:numel(aid)
    site_idx = SID==ii;
    valid_idx = idx_value & site_idx;
    vals = nan(1,numel(midpoints));
    for iv=1:numel(vals)
        d_idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width;
        idx = d_idx & valid_idx;
        vals(iv)=mean(FIRING_RATE(idx));
        
    end
    allSitesRate(ii,:)=vals;
end
figure
subplot(1,2,1)
errF = nanstd(allSitesRate)./sqrt(sum(~isnan(allSitesRate)));
errorbar(nanmean(allSitesRate(:,enough_idx)),midpoints(:,enough_idx),[],[],errF(:,enough_idx),errF(:,enough_idx),'r')
ylabel('Depth')
xlabel('Firing Rate')
subplot(1,2,2)
errC = nanstd(AVG_CORR)./sqrt(sum(~isnan(AVG_CORR)));
errorbar(nanmean(AVG_CORR),corr_midpoints,[],[],errC,errC)
ylim([-2500 1000])
ylabel('Depth')
xlabel('Theta correlation')
        
%% normalize for pxx peaks
aid = unique(SID);
size_idx=~isnan(PXXPeak) & FIRING_RATE>1;
%idx_value = all(PXXSort(:,:,1)>PXXSort(:,:,3),2) & FIRING_RATE'>FR_THRESHOLD;
%size_idx = idx_value';
frac =.1;
peaks_norm = [];
peaks_raw = [];
DEPTH_THIS = [];
col=[];
cntr = 0;

for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    valid_idx_factor = size_idx & idx & DEPTH>-500;
    if nnz(valid_idx)/nnz(idx)<median(frac_dist_coding)
        continue
        
    end
    cntr = cntr+1;
    col = cat(1,col,cntr*ones(nnz(valid_idx),1));
    depth_this = DEPTH(valid_idx);
    peaks_this = PXXPeak(valid_idx);
    %[~,sid]=sort(depth_this,'descend');
    %n=round(numel(sid)*frac);
    %val = median(peaks_this(sid(1:n)));
    peaks_norm_factor = PXXPeak(valid_idx_factor);
    [~,sid]=sort(peaks_norm_factor);
    
    %[~,sid]=sort(peaks_this);
    n=round(frac*numel(sid));
    %val = median(peaks_this(sid(1:n)));
    val = median(peaks_norm_factor(sid(1:n)));
    
    
    peaks_raw = cat(2,peaks_raw,peaks_this);
    peaks_norm=cat(2,peaks_norm,peaks_this/val);
    DEPTH_THIS = cat(2,DEPTH_THIS,depth_this);
%         subplot(1,2,1)
%         histogram(peaks_this,30)
%         subplot(1,2,2)
%         histogram(peaks_this/val,30);
%         pause
%         clf
end
figure
subplot(2,2,1)
scatter(peaks_norm,DEPTH_THIS,30,col,'.')
xlim([0.5 5])
colormap(lines)
subplot(2,2,3)
%scatter(PXXPeak(FIRING_RATE>FR_THRESHOLD),DEPTH(FIRING_RATE>FR_THRESHOLD),30,[.9 .9 .9],'.')
%hold on
scatter(peaks_norm,DEPTH_THIS,30,'.')

xlim([0.5 5])

subplot(2,2,2)
scatter(peaks_raw,DEPTH_THIS,30,col,'.')

%% normalize for xcorr peaks
aid = unique(SID);
size_idx=PEAKValue>QUANT & FIRING_RATE>1;
frac =.1;
peaks_norm = [];
peaks_raw = [];
DEPTH_THIS = [];
col=[];
cntr = 0;
figure
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    valid_idx_factor = size_idx & idx & DEPTH>-500;
    if nnz(valid_idx)/nnz(idx)<median(frac_dist_coding_xcorr)
        continue
        
    end
    cntr = cntr+1;
    col = cat(1,col,cntr*ones(nnz(valid_idx),1));
    depth_this = DEPTH(valid_idx);
    peaks_this = PEAKS(valid_idx);
    %[~,sid]=sort(depth_this,'descend');
    %n=round(numel(sid)*frac);
    %val = median(peaks_this(sid(1:n)));
    peaks_norm_factor = PEAKS(valid_idx_factor);
    [~,sid]=sort(peaks_norm_factor);
    
    %[~,sid]=sort(peaks_this);
    n=round(frac*numel(sid));
    %val = median(peaks_this(sid(1:n)));
    val = median(peaks_norm_factor(sid(1:n)));
    
    
    peaks_raw = cat(2,peaks_raw,peaks_this);
    peaks_norm=cat(2,peaks_norm,peaks_this/val);
    DEPTH_THIS = cat(2,DEPTH_THIS,depth_this);
%         subplot(1,2,1)
%         histogram(peaks_this,30)
%         subplot(1,2,2)
%         histogram(peaks_this/val,30);
%         pause
%         clf
end
figure
subplot(2,2,1)
scatter(peaks_norm,DEPTH_THIS,30,col,'.')
xlim([0.5 5])
colormap(lines)
subplot(2,2,3)
%scatter(PXXPeak(FIRING_RATE>FR_THRESHOLD),DEPTH(FIRING_RATE>FR_THRESHOLD),30,[.9 .9 .9],'.')
%hold on
scatter(peaks_norm,DEPTH_THIS,30,'.')

xlim([0.5 5])

subplot(2,2,2)
scatter(peaks_raw,DEPTH_THIS,30,col,'.')

% %colormap hsv
% figure
% histogram(peaks_norm,150)
% figure
% YBIN = linspace(-2450,1550,21);
% XBIN = linspace(0.4,5,51);
% h = histogram2(peaks_norm,DEPTH_THIS,'YBinEdges',YBIN,'XBinEdges',XBIN,'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
%%
figHandles = findobj('Type', 'figure');
for iF=1:numel(figHandles)
    set(figHandles(iF),'Renderer','Painters')
    saveas(figHandles(iF),sprintf('F:/temp/figures/fig1_baseline_xcorr_%d.pdf',iF))
end