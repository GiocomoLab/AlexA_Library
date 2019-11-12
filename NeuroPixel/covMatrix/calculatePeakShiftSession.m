function [PEAKS,SHIFTS]=calculatePeakShiftSession(data,trials,chunksize,stride,region,stability)
%extract spatial maps

%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
edges=[0:2:400];
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
filt = gausswin(10);
filt = filt/sum(filt);
filt = reshape(filt,[1, numel(filt),1]);

sPF=convn(spatialMap,filt,'same');
% get template trials
template_trials = [4:10];
template = nanmean(sPF(:,:,template_trials),3);

%%
nBins = size(spatialMap,2);
startVec = 1:stride:(nBins-chunksize);
nReps = numel(startVec);
maxlag = 10;
nTrials = size(spatialMap,3);
SHIFTS = nan(nTrials,nReps);
PEAKS = SHIFTS;
repidx = 0;
subset = calc_xcorr_snippet(spatialMap(:,:,template_trials),template,1,200,20);
peaks = max(subset,[],3);
stable_cells = all(peaks>stability,2);
if nnz(stable_cells)<20
    return
end
for iStart = 1:stride:(nBins-chunksize)
    repidx = repidx+1;
    startbin = iStart;
    stopbin = iStart+chunksize;
    
    [xcorrs,lags] = calc_xcorr_snippet(spatialMap,template,startbin,stopbin,maxlag);
    m_xcorr = squeeze(nanmean(xcorrs(stable_cells,:,:),1));
    [peaks,iidx]=max(m_xcorr,[],2);
    shifts = lags(iidx);
    PEAKS(:,repidx)=peaks;
    SHIFTS(:,repidx)=shifts;
end
%%
   

%%
%cm = winter(24);
%figure;hold on; for ii=1:1:24; plot(lags,squeeze(nanmean(xcorrs(:,ii,:),1)),'Color',cm(ii,:)); end

