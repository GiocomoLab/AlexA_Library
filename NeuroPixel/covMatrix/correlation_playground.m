%load data

data = load('F:\NP_DATA\npF4_1025_gaincontrast_2.mat');
%% calculate spatial firing rate map around onset trials (8 pre, 4 post)
region = 'MEC';
onset_trial = 21;
range =-8:3;
trials = onset_trial+range;
binsize=2;
spatialMap=[];
dwell_time=[];
edges=[0:binsize:400];
edges(1)=-.01;
data.posx(data.posx<0)=0;
data.posx(data.posx>=400)=399.00;
for iT=1:length(trials)
    idxVR=data.trial==trials(iT);
    t_time=data.post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=data.sp.st<stop & data.sp.st>=start;
    [spM, dT]=getSpikeMatPosition(data.sp.st(idxNP),data.sp.clu(idxNP),data.posx(idxVR),data.post(idxVR),'edges',edges,'max_clust',max(data.sp.clu)+1);
    spatialMap=cat(3,spatialMap,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
reg = strcmp(data.anatomy.cluster_parent,region);
if iscolumn(reg)
    reg = reg';
end
spatialMap=spatialMap(data.sp.cids+1,:,:);
spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);
%spatialMap = spatialMap(this_stab>0.0 & this_MEC,:,:);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
spatialMap(isnan(spatialMap))=0;
% do spatial smoothing
% filt = gausswin(11);
% filt = filt/sum(filt);
% filt = reshape(filt,[1, numel(filt),1]);
smoothSigma = 4/binsize;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
sPF = repmat(spatialMap,[1,3,1]);
sPF=convn(sPF,filt,'same');
iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
sPF = sPF(:,iidx,:);

spatialMap = sPF;


%calculate Xspatial


%calculate XTXspatial



% for 4 baseline trials, calculate mean response, and then pairwise
% correlation, do some ordering
%%
firing_rates=squeeze(mean(spatialMap(:,:,4:8),3));
firing_rates = normr(firing_rates);
Y=pdist(firing_rates);
Z=linkage(Y);
order = optimalleaforder(Z,Y);
cmat = corr(firing_rates');
figure
subplot(1,2,1)
imagesc(cmat(order,order))

firing_rates=squeeze(mean(spatialMap(:,:,9:end),3));
firing_rates = normr(firing_rates);
cmat = corr(firing_rates');
subplot(1,2,2)
imagesc(cmat(order,order))


%do the same for gain change trials and order with baseline clustering