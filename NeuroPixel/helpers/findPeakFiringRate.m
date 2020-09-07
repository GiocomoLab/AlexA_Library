function [data_out] = findPeakFiringRate(data,ops)

%prepare variables
if isfield(data.anatomy,'parent_shifted')
    reg = data.anatomy.parent_shifted;
else
    reg = data.anatomy.cluster_parent;
end


if isfield(data.anatomy,'parent_shifted')
    sub_reg = data.anatomy.region_shifted;
elseif isfield(data.anatomy,'cluster_region')
    
    sub_reg = data.anatomy.cluster_region;
else
    sub_reg = {};
end

if iscolumn(reg)
    reg = reg';
    sub_reg = sub_reg';
end

if isfield(data.anatomy,'depth_shifted')
    depth = data.anatomy.depth_shifted;
elseif isfield(data.anatomy,'depth')
    depth = data.anatomy.depth;
end
if isfield(data.anatomy,'z2') % for MEC cases
    depth = data.anatomy.tip_distance' - data.anatomy.z2;
end

depth = depth(data.sp.cgs==2);
reg = reg(data.sp.cgs==2);
if ~isempty(sub_reg)
    sub_reg = sub_reg(data.sp.cgs==2);

trials = ops.trials;


%create trial map that only contains numbers for trials to be included
trialMap = nan(1,numel(data.trial_gain));

cntr = 1;
for iT =1:numel(data.trial_gain)
    if ismember(iT,trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end




%recreate data.trial map that resets numbers of included trials and sets
%all else to nan
trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end



[speed,speed_raw]=calcSpeed(data.posx,ops);
if ~isempty(ops.filter)
    speed_raw = conv(speed_raw,ops.filter,'same');
end

spatialMap=[];
dwell_time=[];
edges=ops.edges;
trials = ops.trials;
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


dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
%spatialMap(isnan(spatialMap))=0;
spatialMap = fillmissing(spatialMap,'pchip',2);
smoothSigma = ops.smoothSigma/ops.BinWidth;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
sPF = repmat(spatialMap,[1,3,1]);
sPF=convn(sPF,filt,'same');
iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
sPF = sPF(:,iidx,:);
spMap = sPF;
nT=numel(ops.trials);
speedMat = nan(nT,ops.nBins);
        
        
        discrete_pos = discretize(data.posx,ops.edges);
        cntr=0;
        for iT=ops.trials
            cntr=cntr+1;
            for iB=1:ops.nBins
                idx = data.trial==iT & discrete_pos == iB;
                if nnz(idx)>0
                    val = mean(speed(idx));
                    speedMat(cntr,iB)=val;
                end
            end
        end
        speedMat = fillmissing(speedMat,'pchip',2);






maxBinInd=nan(1,nnz(data.sp.cgs==2));
maxLoc = maxBinInd;
for iC = 1:size(spMap,1)
   
    
    spMap_this = squeeze(spMap(iC,:,:));
    
    tmp=sum(spMap_this,1);
    if nnz(tmp==0)>=2 %more than 2 trials without spikes
        continue
    end
    
    avFR = mean(spMap_this,2);
    [ma,mi]=max(avFR(ops.search_range));
    [ma2,mi2]=max(avFR([ min(ops.search_range)-1 max(ops.search_range)+1] ));
    mi= mi+min(ops.search_range)-1;
    if ma>ma2
        maxLoc(iC)=ops.midpoints(mi);
        maxBinInd(iC)=(mi);
        
    end
end

data_out.region = reg;
data_out.subregion = sub_reg;


data_out.maxBin_index = maxBinInd;
data_out.maxLoc = maxLoc;
data_out.speedMat = speedMat;

data_out.depth = depth;
end



