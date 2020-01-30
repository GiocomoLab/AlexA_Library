function [PEAKS,SHIFTS,speed_mat,stability]=calculatePeakShift(data,trials,ops)

spatialMap=[];
dwell_time=[];
edges=[0:ops.binsize:400];
edges(1)=-.01;
[speed,speed_mat] = calc_speed(data.posx,data.trial,trials,edges);
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

spatialMap=spatialMap(data.sp.cids+1,:,:);


dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
spatialMap(isnan(spatialMap))=0;
smoothSigma = ops.smoothSigma/ops.binsize;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
sPF = repmat(spatialMap,[1,3,1]);
sPF=convn(sPF,filt,'same');
iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
sPF = sPF(:,iidx,:);
% get template trials

spatialMap = sPF;


fr_mat = squeeze(nanmean(spatialMap));
spatialMap=spatialMap(data.sp.cgs==2,:,:);


template1_trials = ops.template_trials;
template1={};
cntr = 0;



for ii=template1_trials
    cntr = cntr+1;
    tmp_idx =setdiff(template1_trials,ii);

    tmp = nanmean(spatialMap(:,:,tmp_idx),3);
    template1{cntr}=tmp;
end
template1{numel(template1_trials)+1}=nanmean(spatialMap(:,:,template1_trials),3);


%%
nBins = size(spatialMap,2);
startVec = ops.stride_start:ops.stride:(nBins-ops.chunksize+1);
nReps = numel(startVec);
maxlag = 10;
nTrials = size(spatialMap,3);
SHIFTS = nan(nTrials,nReps);
PEAKS = SHIFTS;
repidx = 0;
%stable cells: have a peak xcorr across the whole thing of greater than
%thresh for each trial
idx = find(triu(true(numel(ops.template_trials)),1));
stability = zeros(1,size(spatialMap,1));
for ii=1:numel(stability)
    tmp = squeeze(spatialMap(ii,:,ops.template_trials));
    tmp = corr(tmp);
    stability(ii)=mean(tmp(idx));
end



trial2templateMap=zeros(1,size(spatialMap,3));
cntr=1;
aa=numel(template1_trials)+1;
for iT = 1:numel(trial2templateMap)
    if ~ismember(iT,ops.template_trials)
        trial2templateMap(iT)=aa;
    else
        trial2templateMap(iT)=cntr;
        cntr=cntr+1;
    end
end
nCells = size(spatialMap,1);

SHIFTS=zeros(nCells,nTrials,nReps);
PEAKS = SHIFTS;
repidx = 0;
for iStart = ops.stride_start:ops.stride:(nBins-ops.chunksize)
    repidx = repidx+1;
    startbin = iStart;
    stopbin = iStart+ops.chunksize;

    [xcorrs,lags] = calc_xcorr_snippet(spatialMap,template1,startbin,stopbin,maxlag,trial2templateMap);
    [ma,mi]=max(xcorrs,[],3);
    PEAKS(:,:,repidx)=ma;
    SHIFTS(:,:,repidx)=lags(mi);
end
end

function [speed,speed_mat] = calc_speed(posx,trial_pos,trials,edges)
speed = diff(posx)/0.02;

% throw out extreme values and interpolate
speed(speed > 150) = NaN;
speed(speed<-5) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = [0;speed];
posx(posx<0)=0;
posx(posx>=400)=399.99;

d_loc = discretize(posx,edges);
speed_mat = nan(numel(trials,numel(edges)-2));
for iT=1:numel(trials)
    idxT = trial_pos == trials(iT);

    for iB=1:numel(edges)-1
        idx= idxT & d_loc==iB;
        speed_mat(iT,iB)=mean(speed(idx));
    end
end

end