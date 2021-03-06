function [PEAKS,SHIFTS,XTX,n_units,YYT,speed_mat,fr_mat,MEAN_XCORRS]=calculatePeakShiftSession(data,trials,ops)
%extract spatial maps

%trials = trials(trial_gain == 1 & trial_contrast == 100);
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
%cellIDX=find(sp.cgs>=1);
if isfield(data.anatomy,'parent_shifted')
    reg = startsWith(data.anatomy.parent_shifted,ops.region);
else
reg = startsWith(data.anatomy.cluster_parent,ops.region);
end
if iscolumn(reg)
    reg = reg';
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


%pp=bsxfun(@rdivide,spatialMap,mm);
fr_mat = squeeze(nanmean(spatialMap));
spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);


template1_trials = ops.template_trials;
template1={};
cntr = 0;

%depth = data.anatomy.tip_distance(data.sp.cgs==2 & reg);
%sel_idx = depth<median(depth);
%spatialMap = spatialMap(sel_idx,:,:);

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
% fit a glm to the data
% clu_list = data.sp.cids(data.sp.cgs==2 & reg);
% tr = max(min(trials)-4,2):min(trials)+5;
% glmData = fitGLM(data,tr,clu_list);
% glmscore = nan(2,numel(glmData));
% for iC=1:numel(glmData)
%     if iscell(glmData(iC).allModelTestFits)
%     glmscore(:,iC)=mean(glmData(iC).allModelTestFits{1});
%     end
% end
stable_cells = stability>ops.stability_threshold;
%stable_cells = glmscore(2,:)>0 & glmscore(2,:)>glmscore(1,:);

% subset = calc_xcorr_snippet(spatialMap(:,:,1:10),template1,1,200,20);
% peaks = max(subset,[],3);
% stable_cells = all(peaks>stability,2);
XTX = nan(nTrials*nBins);
YYT = nan(nTrials,nTrials);
n_units = [nnz(stable_cells),numel(stable_cells)];
MEAN_XCORRS = nan(nTrials,2*maxlag+1,nReps);
if nnz(stable_cells)<=.2*numel(stable_cells)
    return
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

%spatialMap = spatialMap-mean(spatialMap,2);
for iStart = ops.stride_start:ops.stride:(nBins-ops.chunksize)
    repidx = repidx+1;
    startbin = iStart;
    stopbin = iStart+ops.chunksize;
    %xcorr relative to template trials: 4 before gain onset for all except
    %these 4
    [xcorrs,lags] = calc_xcorr_snippet(spatialMap,template1,startbin,stopbin,maxlag,trial2templateMap);
    %[xcorrs2,~] = calc_xcorr_snippet(spatialMap,template2,startbin,stopbin,maxlag);
    
    m_xcorr = squeeze(nanmean(xcorrs(stable_cells,:,:),1));
    MEAN_XCORRS(:,:,repidx)=m_xcorr;
    [peaks,iidx]=max(m_xcorr,[],2);
    shifts = lags(iidx);

    %[pp,pps] = max(xcorrs,[],3);
    %peaks=nanmean(pp);
    %midx = round(nanmean(pps));
    %shifts = lags(midx);
    
    %m_xcorr2 = squeeze(nanmean(xcorrs2(stable_cells,:,:),1));
    %[peaks2,iidx2]=max(m_xcorr2,[],2);
    %peaks(template1_trials)=peaks2(template1_trials);
    %iidx(template1_trials)=iidx2(template1_trials);
    PEAKS(:,repidx)=peaks;
    SHIFTS(:,repidx)=shifts;
end
% trial by trial similarity, all cells
% Y=zeros(size(spatialMap,3),nnz(stable_cells)*size(spatialMap,2));
% for iT=1:size(Y,1)
%     tmp = squeeze(spatialMap(stable_cells,:,iT))';
%     tmp = reshape(tmp,1,[]);
%     Y(iT,:)=tmp;
% end
% 
% Y=Y-mean(Y,1);
% Y=normr(Y);
% 
% YYT = Y*Y';
spatialMap = spatialMap(stable_cells,:,:);
correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
for iC=1:size(spatialMap,1)
    tmp=corr(squeeze(spatialMap(iC,:,:)));
    correlation_All(:,:,iC)=tmp;
end
YYT = squeeze(nanmean(correlation_All,3));
% position by position covariance matrix, stable cells
%spatialMap = spatialMap(stable_cells,:,:);

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



%%
%cm = winter(24);
%figure;hold on; for ii=1:1:24; plot(lags,squeeze(nanmean(xcorrs(:,ii,:),1)),'Color',cm(ii,:)); end

