function [XTX,trial,gain,contrast]=getCovMatrix(data,cell_info_session,trials)

local_stab = cell_info_session.LocalStability;
numblocks = floor(32/4);
local_stab = mean(local_stab(:,1:numblocks),2);
this_stab = local_stab;
this_reg = cell_info_session.BrainRegion;
this_MEC = strcmp(this_reg,'MEC');



%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
edges=[0:5:400];
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
spatialMap=spatialMap(data.sp.cids+1,:,:);
spatialMap=spatialMap(data.sp.cgs==2,:,:);
spatialMap = spatialMap(this_stab>0.0 & this_MEC,:,:);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
spatialMap(isnan(spatialMap))=0;
%%

X=zeros(size(spatialMap,1),size(spatialMap,2)*size(spatialMap,3));
filt = gausswin(15);
filt = filt/sum(filt);
for iC = 1:size(spatialMap,1)
    tmp = spatialMap(iC,:,:);
    
    X(iC,:)=tmp(:);
end
X(isnan(X))=0;
X=conv2(X,filt','same');
X=X-mean(X,2);
X=normc(X);

XTX = X'*X;
trial = trials;
gain = data.trial_gain(trials);
contrast = data.trial_contrast(trials);
end
