%load data

data = load('/Users/attialex/Desktop/data/npF3_1018_contrasttrack_gainchanges_1.mat');
%% calculate spatial firing rate map around onset trials (8 pre, 4 post)
region = 'MEC';
onset_trial = 21;
range =-8:3;
trials = onset_trial+range;
binsize=2;
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
reg = strcmp(data.anatomy.cluster_parent,region);
if iscolumn(reg)
    reg = reg';
end
spatialMap=spatialMap(data.sp.cids+1,:,:);
spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);
good_cells = data.sp.cids(data.sp.cgs==2 & reg);
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

spatialMap = sPF;


%calculate Xspatial


%calculate XTXspatial



% for 4 baseline trials, calculate mean response, and then pairwise
% correlation, do some ordering
%%
normc = @(m) sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
normr = @(m) sqrt(m.^2 ./ sum(m.^2,2)) .* sign(m);
%normr = @(m) sqrt(ones./(sum((m.*m)')))'*ones(1,mc).*m;
%% 
X=zeros(size(spatialMap,1),size(spatialMap,2)*8);

for iC = 1:size(X,1)
    tmp = spatialMap(iC,:,5:12);
    
    X(iC,:)=tmp(:);
end
X(isnan(X))=0;
X=X-mean(X,2);
X=normc(X);

XTX = X'*X;
figure
imagesc(XTX)

%% correlation of spatial firing maps
spatial_firing_rates_pre=squeeze(mean(spatialMap(:,:,4:8),3));
spatial_firing_rates_pre = normr(spatial_firing_rates_pre);
b=~isnan(spatial_firing_rates_pre(:,1));
spatial_firing_rates_pre = spatial_firing_rates_pre(b,:);
Y=pdist(spatial_firing_rates_pre);
Z=linkage(Y);
order = optimalleaforder(Z,Y);
cmat = corr(spatial_firing_rates_pre');
figure
subplot(1,2,1)
imagesc(cmat(order,order))

spatial_firing_rates_post=squeeze(mean(spatialMap(:,:,9:end),3));
spatial_firing_rates_post = normr(spatial_firing_rates_post);
spatial_firing_rates_post=spatial_firing_rates_post(b,:);
cmat_post = corr(spatial_firing_rates_post');
subplot(1,2,2)
imagesc(cmat_post(order,order))

%% now calculate temporal correlations

time_bin = 0.1;
time_range=data.post(ismember(data.trial,[21:24]));
start_time = min(time_range);
stop_time = max(time_range);
duration = stop_time-start_time;
start_time = start_time-duration;
time_vec = start_time:time_bin:stop_time;
middle = min(data.post(ismember(data.trial,21)));
[~,middle_bin] = min(abs(middle-time_vec));

X = zeros(numel(good_cells),numel(time_vec)-1);

for iC=1:numel(good_cells)
    X(iC,:)=histcounts(data.sp.st(data.sp.clu==good_cells(iC)),time_vec);
end

filt=gausswin(1)';
filt = filt/sum(filt);
X=conv2(X,filt,'same');
X=X-mean(X,2);
%X=normc(X);
X=normr(X);
b=~isnan(X(:,1));
X=X(b,:);
%%
Y = pdist(X(:,1:middle_bin));
Z = linkage(Y);
order = optimalleaforder(Z,Y);
figure
tmp_pre = corr(X(:,1:middle_bin)');
subplot(1,2,1)
imagesc(tmp_pre(order,order))

tmp_post = corr(X(:,middle_bin:end)');
subplot(1,2,2)
imagesc(tmp_post(order,order))
%%
cc = tmp_pre-tmp_post;
%extract upper half
diff_cc = cc(triu(true(size(cc)),1));
