function [spatial_similarity,spatial_similarity_smooth] = compute_spatial_similarity(dataset)
sp = dataset.sp;
trial = dataset.trial;
post = dataset.post;
posx = dataset.posx;
trials=[1:max(trial)];

spatialMap=[];
dwell_time=[];
edges=[0:5:405];
edges(1)=-.01;
posx(posx<0)=0;
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.clu)+1);
    spatialMap=cat(3,spatialMap,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
spatialMap=spatialMap(sp.cids+1,:,:);
spatialMap=spatialMap(sp.cgs==2,:,:);
spatialMap=spatialMap(:,1:end-1,:);
dwell_time=dwell_time(:,1:end-1);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
%fill in nan values
if nnz(isnan(spatialMap))>0
    spatialMap = fillmissing(spatialMap,'linear',2);
end
%smooth
 win = gausswin(11);
 gauss_filter = win/sum(win);
 gauss_filter = reshape(gauss_filter,[1,numel(gauss_filter),]);
 spatialMap_smooth = convn(spatialMap,gauss_filter,'same');

correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
correlation_smooth = correlation_All;
for iC=1:size(spatialMap,1)
    tmp=corr(squeeze(spatialMap(iC,:,:)));
    correlation_All(:,:,iC)=tmp;
    tmp=corr(squeeze(spatialMap_smooth(iC,:,:)));
    correlation_smooth(:,:,iC)=tmp;
end
spatial_similarity=correlation_All;
spatial_similarity_smooth = correlation_smooth;
end
