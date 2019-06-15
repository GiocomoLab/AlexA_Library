function [spatialMap,spatialMap_smoothed]=get_spatial_map(dataset,trialset)
sp = dataset.sp;
trial = dataset.trial;
post = dataset.post;
posx = dataset.posx;
posx(posx>400)=400;
if nargin ==1
trialset=[1:max(trial)];
end

spatialMap=[];
dwell_time=[];
edges=[0:4:400];
edges(1)=-.01;
posx(posx<0)=0;
for iT=1:length(trialset)
    idxVR=trial==trialset(iT);
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
%spatialMap=spatialMap(:,1:end-1,:);
%dwell_time=dwell_time(:,1:end-1);
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
 win = gausswin(9);
 gauss_filter = win/sum(win);
 gauss_filter = reshape(gauss_filter,[1,numel(gauss_filter),]);
 spatialMap_smoothed = convn(spatialMap,gauss_filter,'same');
 
 %tuning_curve = [];
 %tuning_curve_smoothed = [];
end