function [meanSpatialMap,meanSpatialMapNormalized,nunits,respIncreasing] = getSpatialMap(data,region)

trials = find(data.trial_gain == 1 & data.trial_contrast == 100);

spatialMap = calculateSpatialFiringMap(data,trials,2);

if isfield(data.anatomy,'parent_shifted')
    reg = startsWith(data.anatomy.parent_shifted,region);
else
reg = startsWith(data.anatomy.cluster_parent,region);
end
if iscolumn(reg)
    reg = reg';
end

spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);

meanSpatialMap = squeeze(nanmean(spatialMap,3));
meanSpatialMap = nanmean(meanSpatialMap,1);

meanE = squeeze(nanmean(spatialMap(:,:,2:2:end),3));
meanO = squeeze(nanmean(spatialMap(:,:,1:2:end),3));

respLate=squeeze(nanmean(spatialMap(:,141:190,2:2:end),2));
respEarly = squeeze(nanmean(spatialMap(:,1:50,2:2:end),2));
Pvals = nan(size(spatialMap,1),1);
for iC=1:numel(Pvals)
    x=squeeze(respLate(iC,:)-respEarly(iC,:));
    tmp = signrank(x,0,'tail','right');
    Pvals(iC)=tmp;
end
respIncreasing = nanmean(meanO(Pvals<0.05,:));
X=zeros(size(spatialMap,1),size(spatialMap,2)*size(spatialMap,3));

for iC = 1:size(X,1)
    tmp = spatialMap(iC,:,:);
    
    X(iC,:)=tmp(:);
end
X(isnan(X))=0;
%X=X-mean(X,2);
%X=normc(X);
X=zscore(X,0,2);

spN = spatialMap;
nunits = size(spatialMap,1);
for iC = 1:size(X,1)
    
    for iT = 1:size(spatialMap,3)
        idx = ((iT-1)*200+1):iT*200;
        tmp = X(iC,idx);
        spN(iC,:,iT)=tmp;
    end
end
meanSpatialMapNormalized = squeeze(nanmean(spN,3));
meanSpatialMapNormalized = nanmean(meanSpatialMapNormalized,1);



end

function spatialMap = calculateSpatialFiringMap(data,trials,binsize)

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

spatialMap=spatialMap(data.sp.cids+1,:,:);
%spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);
%spatialMap = spatialMap(this_stab>0.0 & this_MEC,:,:);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
spatialMap(isnan(spatialMap))=0;
end