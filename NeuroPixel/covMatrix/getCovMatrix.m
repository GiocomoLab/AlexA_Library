function [XTX,YYT,trial,gain,contrast,nunits]=getCovMatrix(data,region,trials,binsize,template_trials,stability_threshold)
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
filt = gauss_filter;
%filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
% sPF = repmat(spatialMap,[1,3,1]);
% sPF=convn(sPF,filt,'same');
% iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
% sPF = sPF(:,iidx,:);
%get template trials

%spatialMap = sPF;
nunits = size(spatialMap,1);
%%

idx = find(triu(true(numel(template_trials)),1));
stability = zeros(1,size(spatialMap,1));
for ii=1:numel(stability)
    tmp = squeeze(spatialMap(ii,:,template_trials));
    tmp = corr(tmp);
    stability(ii)=mean(tmp(idx));
end
    
stable_cells = stability>stability_threshold;

X=zeros(size(spatialMap,1),size(spatialMap,2)*size(spatialMap,3));
for iC = 1:size(spatialMap,1)
    tmp = spatialMap(iC,:,:);
    
    X(iC,:)=tmp(:);
end
X(isnan(X))=0;
X=conv2(X,filt','same');
X=X(stable_cells,:);
X=X-mean(X,2);
X=normc(X);

XTX = X'*X;
trial = trials;
gain = data.trial_gain(trials);
contrast = data.trial_contrast(trials);


filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
sPF = repmat(spatialMap,[1,3,1]);
sPF=convn(sPF,filt,'same');
iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
sPF = sPF(:,iidx,:);

spatialMap = sPF;

Y=zeros(size(spatialMap,3),size(spatialMap,1)*size(spatialMap,2));
for iT=1:size(Y,1)
    tmp = squeeze(spatialMap(:,:,iT))';
    tmp = reshape(tmp,1,[]);
    Y(iT,:)=tmp;
end

Y=Y-mean(Y,1);
Y=normr(Y);

YYT = Y*Y';

% correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
% for iC=1:size(spatialMap,1)
%     tmp=corr(squeeze(spatialMap(iC,:,:)));
%     correlation_All(:,:,iC)=tmp;
%     diagAll(iC,:)=diag(tmp,1);
% end

end
