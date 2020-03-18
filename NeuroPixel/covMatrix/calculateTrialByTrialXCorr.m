function [corrMat,shiftMat,stability,spatialMap]=calculateTrialByTrialXCorr(data,ops,bins2correlate)
%extract spatial maps
if nargin ==2
    bins2correlate = 1:ops.nBins;
end
nCells = nnz(data.sp.cgs==2);
if isvector(bins2correlate)
    if iscolumn(bins2correlate)
        bins2correlate=bins2correlate';
    end
    bins2correlate = repmat(bins2correlate,nCells,1);
end
%trials = trials(trial_gain == 1 & trial_contrast == 100);
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
% get template trials

spatialMap = sPF;


idx = (triu(true(numel(ops.trials)),1));
stability = nan(1,size(spatialMap,1));
for ii=1:numel(stability)
    if any(isnan(bins2correlate(ii,:)))
        continue
        
    end
    tmp = squeeze(spatialMap(ii,bins2correlate(ii,:),:));
    tmp = corr(tmp);
    stability(ii)=nanmean(tmp(idx));
end




corrMat = nan(numel(stability),numel(trials),numel(trials));
shiftMat = nan(numel(stability),numel(trials),numel(trials));
shift_all = -ops.maxLag:ops.BinWidth:ops.maxLag;
for cellIDX = 1:size(spatialMap,1)
    if any(isnan(bins2correlate(cellIDX,:)))
        continue
        
    end
    
    mS=squeeze(spatialMap(cellIDX,:,:));
    
    mS=mS-mean(mS);
    mS = mS(bins2correlate(cellIDX,:),:);
    %mS=zscore(mS);
    %%
    [xcorr_this,~]=xcorr(mS,'coeff',ops.maxLag/ops.BinWidth);
    [xcorr_this,max_idx] = max(xcorr_this,[],1); % take max over lags
    shift_this = shift_all(max_idx);
    xcorr_this = reshape(xcorr_this,numel(trials),numel(trials));
    shiftMat(cellIDX,:,:) = reshape(shift_this,numel(trials),numel(trials));
    corrMat(cellIDX,:,:) = xcorr_this-diag(diag(xcorr_this)); % subtract diags
    
end




end
%%
