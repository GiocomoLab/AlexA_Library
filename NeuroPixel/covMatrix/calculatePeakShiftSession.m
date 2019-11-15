function [PEAKS,SHIFTS,XTX]=calculatePeakShiftSession(data,trials,chunksize,stride_start,stride,region,stability,binsize)
%extract spatial maps

%trials = trials(trial_gain == 1 & trial_contrast == 100);
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
if isfield(data.anatomy,'parent_shifted')
    reg = startsWith(data.anatomy.parent_shifted,region);
else
reg = startsWith(data.anatomy.cluster_parent,region);
end
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
% get template trials

spatialMap = sPF;
template1_trials = [1:8];
template2_trials = [1:5];
template1={};
cntr = 0;
for ii=template1_trials
    cntr = cntr+1;
    tmp_idx =setdiff(template1_trials,ii);

    tmp = nanmean(spatialMap(:,:,tmp_idx),3);
    template1{cntr}=tmp;
end
template1{numel(template1_trials)+1}=nanmean(spatialMap(:,:,template1_trials),3);
%template1 = nanmean(spatialMap(:,:,template1_trials),3);
%template1 = template1-mean(template1,2);
template2 = nanmean(spatialMap(:,:,template2_trials),3);

%%
nBins = size(spatialMap,2);
startVec = stride_start:stride:(nBins-chunksize+1);
nReps = numel(startVec);
maxlag = 10;
nTrials = size(spatialMap,3);
SHIFTS = nan(nTrials,nReps);
PEAKS = SHIFTS;
repidx = 0;
%stable cells: have a peak xcorr across the whole thing of greater than
%thresh for each trial
idx = find(triu(true(8),1));
stability = zeros(1,size(spatialMap,1));
for ii=1:numel(stability)
    tmp = squeeze(spatialMap(ii,:,2:10));
    tmp = corr(tmp);
    stability(ii)=mean(tmp(idx));
end
    
stable_cells = stability>.2;
% subset = calc_xcorr_snippet(spatialMap(:,:,1:10),template1,1,200,20);
% peaks = max(subset,[],3);
% stable_cells = all(peaks>stability,2);
XTX = nan(nTrials*nBins);
if nnz(stable_cells)<.2*numel(stable_cells)
    return
end
%spatialMap = spatialMap-mean(spatialMap,2);
for iStart = stride_start:stride:(nBins-chunksize)
    repidx = repidx+1;
    startbin = iStart;
    stopbin = iStart+chunksize;
    %xcorr relative to template trials: 4 before gain onset for all except
    %these 4
    [xcorrs,lags] = calc_xcorr_snippet(spatialMap,template1,startbin,stopbin,maxlag);
    %[xcorrs2,~] = calc_xcorr_snippet(spatialMap,template2,startbin,stopbin,maxlag);
    m_xcorr = squeeze(nanmean(xcorrs(stable_cells,:,:),1));
    [peaks,iidx]=max(m_xcorr,[],2);
    %m_xcorr2 = squeeze(nanmean(xcorrs2(stable_cells,:,:),1));
    %[peaks2,iidx2]=max(m_xcorr2,[],2);
    %peaks(template1_trials)=peaks2(template1_trials);
    %iidx(template1_trials)=iidx2(template1_trials);
    shifts = lags(iidx);
    PEAKS(:,repidx)=peaks;
    SHIFTS(:,repidx)=shifts;
end
spatialMap = spatialMap(stable_cells,:,:);
X=zeros(nnz(stable_cells),size(spatialMap,2)*size(spatialMap,3));

for iC = 1:size(X,1)
    tmp = spatialMap(iC,:,:);
    
    X(iC,:)=tmp(:);
end
X(isnan(X))=0;
X=X-mean(X,2);
X=normc(X);

XTX = X'*X;

end
%%
   

%%
%cm = winter(24);
%figure;hold on; for ii=1:1:24; plot(lags,squeeze(nanmean(xcorrs(:,ii,:),1)),'Color',cm(ii,:)); end

