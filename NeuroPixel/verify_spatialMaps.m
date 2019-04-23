addpath(genpath('C:\code\AlexA_Library'))

trials=[1:max(trial)];
%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
edges=[0:2:410];
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
spatialMap=spatialMap(sp.cgs>1,:,:);
spatialMap=spatialMap(:,1:end-1,:);
dwell_time=dwell_time(:,1:end-1);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
%%
spMap=shiftdim(spatialMap,1);
view_stack(spMap)
%%
spatialMap=[];
dwell_time=[];
edges=[0:10:410];
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

correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
diagAll=zeros(size(spatialMap,1),size(spatialMap,3)-1);
for iC=1:size(spatialMap,1)
    tmp=corr(squeeze(spatialMap(iC,:,:)));
    correlation_All(:,:,iC)=tmp;
    diagAll(iC,:)=diag(tmp,1);
end

figure
pp=nanmean(correlation_All,3);
imagesc(nanmean(correlation_All,3),[-.1 0.3])
title('mean over all corr')
transitions=find(diag(pp,1)<.10);
hold on
plot(transitions+1,transitions,'ro')