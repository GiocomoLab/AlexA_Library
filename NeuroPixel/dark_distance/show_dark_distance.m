files = dir('Z:\giocomo\attialex\distance_coding\data/np*dark*');
PEAKS=[];
PEAKValue=[];
QUANT = [];
DEPTH = [];
PXXPeak = [];
AMP=[];
centers = linspace(1000,-2000,30)
ACG_LIST = {numel(centers,1)};
SID=[];
SMALLMOD = zeros(numel(files));
for iL = 1:numel(centers);
    ACG_LIST{iL}=[];
end
for iF=1:numel(files)
    load(fullfile(files(iF).folder,files(iF).name))
    [a,b]=sort(mec_depth);
%     imagesc(ACG(b,:),[0 0.4])
%     [x,y]=ginput()
%     SMALLMOD(iF)=x(end);
    v_idx =1./f_vec > 25 &  1./f_vec<600;

    [ma,mi]=max(PXX(v_idx,:));
    offset = strfind(v_idx,[0 1]);
    mi = mi+offset;
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
%% normalize for xcorr peaks
aid = unique(SID);
size_idx=PEAKValue>(QUANT+.00);
frac =.2;
peaks_norm = [];
figure
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    depth_this = DEPTH(valid_idx);
    peaks_this = PEAKS(valid_idx);
    [~,sid]=sort(depth_this,'descend');
    n=round(numel(sid)*frac);
    val = median(peaks_this(sid(1:n)));
    peaks_norm=cat(2,peaks_norm,peaks_this/val);
%     subplot(1,2,1)
%     histogram(peaks_this)
%     subplot(1,2,2)
%     histogram(peaks_this/val,20);
%     pause
%     clf
end
figure
scatter(peaks_norm,DEPTH(size_idx))
figure
histogram(peaks_norm,100)
figure

h = histogram2(peaks_norm,DEPTH(size_idx),[12 12],'FaceColor','flat');

%% normalize for pxx peaks
aid = unique(SID);
size_idx=~isnan(PXXPeak);
frac =.2;
peaks_norm = [];
figure
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    depth_this = DEPTH(valid_idx);
    peaks_this = PXXPeak(valid_idx);
    %[~,sid]=sort(depth_this,'descend');
    %n=round(numel(sid)*frac);
    %val = median(peaks_this(sid(1:n)));
    [~,sid]=sort(peaks_this);
    n=round(frac*numel(sid));
    val = mean(peaks_this(sid(1:n)));
    peaks_norm=cat(2,peaks_norm,peaks_this/val);
%     subplot(1,2,1)
%     histogram(peaks_this)
%     subplot(1,2,2)
%     histogram(peaks_this/val,30);
%     pause
%     clf
end
figure
scatter(peaks_norm,DEPTH(size_idx))
figure
histogram(peaks_norm,100)
figure

h = histogram2(peaks_norm,DEPTH(size_idx),[12 12],'FaceColor','flat');
%%
figure
histogram(PEAKS(PEAKValue>(QUANT+.00))*5,50)
%%
figure
idx = PEAKValue>(QUANT);
scatter(PEAKS(idx)*5,DEPTH(idx))

figure
scatter(PXXPeak,DEPTH,'.')
figure
histogram(PXXPeak,100)

%%
% figure
% valididx = ~isnan(PXXPeak) & ~isnan(DEPTH);
% scatplot(PXXPeak(valididx),DEPTH(valididx))

%%
figure
subplot(1,2,1)
idx = PEAKValue>(QUANT+.0) & AMP>0.0;
h = histogram2(PEAKS(idx)*5,DEPTH(idx),[50 50],'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
set(gca,'CLim',[0 15])
xlabel('Peak location [cm]')
ylabel('Distance from MEC entry [um]')
colorbar
set(gcf,'Render','Painters')
subplot(1,2,2)

aid = unique(SID);
frac=zeros(1,numel(aid));
for ii=1:numel(aid)
    idx = SID==ii;
    valid_idx = size_idx & idx;
    frac(ii)=sum(valid_idx)/sum(idx);
end
plotSpread(frac','SpreadWidth',0.2)

boxplot(frac)
%%
figure
idx = PEAKValue>(QUANT+.01);
h = histogram2(PXXPeak,DEPTH,[100 100],'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
set(gca,'CLim',[0 6])
figure
idx = PEAKValue>(QUANT+.01);
h = histogram2(PXXPeak,DEPTH,50:10:300,linspace(-2500,1000,100),'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
set(gca,'CLim',[1 10])
%%
p=zeros(12,101);
for iC=1:12
    p(iC,:)=mean(ACG_LIST{iC});
end
figure
imagesc(p,[0 0.3])
