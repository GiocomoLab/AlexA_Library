function data_out = calculateTrialByTrialXCorr(data,trials)

%housekeeping
binsize=2;
good_cells = data.sp.cids(data.sp.cgs==2);
data_out = struct();
if ~isfield(data,'anatomy')
    return
end
%calculate spatial maps

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
    reg = data.anatomy.parent_shifted;
else
    reg = data.anatomy.cluster_parent;
end
if iscolumn(reg)
    reg = reg';
end
spatialMap=spatialMap(data.sp.cids+1,:,:);
spatialMap = spatialMap(data.sp.cgs==2,:,:);
reg = reg(data.sp.cgs==2);

dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
spatialMap(isnan(spatialMap))=0;
smoothSigma = 4/binsize;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
sPF = repmat(spatialMap,[1,3,1]);
sPF=convn(sPF,filt,'same');
iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
sPF = sPF(:,iidx,:);

spatialMap = sPF;

% calculate stability
idx = (triu(true(numel(trials)),1));
stability = zeros(1,size(spatialMap,1));
for ii=1:numel(stability)
    tmp = squeeze(spatialMap(ii,:,:));
    tmp = corr(tmp);
    stability(ii)=mean(tmp(idx));
end

data_out.stability=stability;
data_out.region = reg;

% calculate speed
posWindow= [300 390];
%posWindow= [120 200];
[speed,~] = calc_speed(data.posx,data.trial,trials,edges);

posBin = [find(edges == posWindow(1)), find(edges == posWindow(2))];
trial_speed = zeros(1,max(data.trial));
for ii = 1:max(data.trial)
    idx = data.posx>posWindow(1) & data.posx<posWindow(2) & data.trial == ii;
    tmp = mean(speed(idx));
    
    trial_speed(ii)=tmp;
end


% calculate pairwise xcorrs for all cells
nC=size(spatialMap,1);
nPairs = length(trials)*(length(trials)-1)/2;
nDat = 3;
allDelays = zeros(nC,nPairs,nDat);
for cellIDX = 1:size(spatialMap,1)
 
    mS=squeeze(spatialMap(cellIDX,:,:));
    
    mS=mS-mean(mS);
    
    %%
    [rr,lags]=xcorr(mS(posBin(1):posBin(2),:),'coeff',10);
    nR=size(mS,2);
    delay_list=zeros(nR*(nR-1)/2,3); %delay,value, speed_diff, 
    idx = 0;
    nT = numel(trials);
    for iT = 1:nT
        for jT=(iT+1):nT
            idx = idx+1;
            c_idx = (iT-1)*nT+jT;
            [mm,midx]=max(rr(:,c_idx));
            if ismember(midx,[1 size(rr,1)])
                midx = nan;
            end
            
            tmp_d = min(lags)-1;
            delay_list(idx,:)=[midx+tmp_d,mm,trial_speed(trials(iT))-trial_speed(trials(jT))];
        end
    end
    allDelays(cellIDX,:,:)=delay_list;
end
data_out.allDelays = allDelays;
data_out.goodCells = data.sp.cids(data.sp.cgs==2);
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
