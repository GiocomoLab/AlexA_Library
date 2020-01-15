load('F:\NP_DATA\npI5_0414_dark_1.mat')
speed = diff(posx);

% throw out extreme values and interpolate

speed(speed<-10) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = [0;speed];
distance = cumsum(speed);
distance_trial = ones(size(post));
%%


binsize=10;
edges=[0:binsize:(max(distance)+binsize)];

idxVR=true(size(distance_trial));
t_time=post(idxVR);
start=min(t_time);
stop=max(t_time);
idxNP=sp.st<stop & sp.st>=start;
[spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),distance,post(idxVR),'edges',edges,'max_clust',max(sp.clu)+1);
%%
spatialMap=bsxfun(@rdivide,spM(sp.cids+1,:,:),dT);
%spatialMap = spM;
spatialMap(isnan(spatialMap))=0;

smoothSigma = 20/binsize;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

%spatialMap = conv2(spatialMap,gauss_filter','same');

%%
good_cells = sp.cgs==2 & strcmp(anatomy.cluster_parent,'MEC')';
depth = anatomy.tip_distance(good_cells);
[~,sid]=sort(depth,'descend');

zSpatialMap = zscore(spatialMap(good_cells,:),0,[2]);
figure
imagesc(zSpatialMap(sid,:),[0 2])
%%
lags = 1000;
nUnits = size(zSpatialMap,1);
ACG = zeros(nUnits,lags+1);
ACG_raw = ACG;
spatialMapr = spatialMap(good_cells,:);
for ii=1:nUnits
    tmp=xcorr(zSpatialMap(ii,:),1000,'coeff');
    ACG(ii,:)=tmp(1001:end);
    tmp = xcorr(spatialMapr(ii,:),1000,'coeff');
    ACG_raw(ii,:)=tmp(1001:end);
end
ACG=ACG;
figure
subplot(1,2,1)
imagesc(ACG(sid,1:50),[0 0.4])
subplot(1,2,2)
imagesc(ACG_raw(sid,1:50))
%%

[pxx,f]=pwelch(zSpatialMap',[],[],[],100/binsize);
[aa,ii]=max(pxx,[],1);
figure
scatter(log(f(ii)),depth)
figure
for ii=1:100
    subplot(1,2,1)
    a=pwelch(zSpatialMap(ii,:)',[],[],[],100/binsize);
    b=pwelch(spatialMapr(ii,:)',[],[],[],100/binsize);
    plot(f,log(a))
    [ma,mi]=max(a);
    
    hold on
    %plot(f,log(b))
    subplot(1,2,2)
    plot(ACG(ii,1:50))
    if ma>upper_bound_pxx(ii,mi)
    xline(1/f(mi)*10,'r')
    end
    hold on
    %plot(ACG_raw(ii,1:50))
    pause
    clf
    
end
%%
PI=[];
figure
imagesc(ACG(sid,1:50),[0 0.4])
hold on
for ii=1:173
    [p,ip]=findpeaks(ACG(sid(ii),1:50),'MinPeakHeight',0,'MinPeakProminence',.0);
    if ~isempty(p)
            PI(end+1)=ip(1);

        if  p(1)>upper_bound(sid(ii),ip(1))
            plot(ip(1),ii,'kx')
        else
            plot(ip(1),ii,'rx')
        end
    end
end

%%
n_spikes = numel(sp.st);
nUnits = size(zSpatialMap,1);
lags=100;
n_it = 100;
abs_min = min(sp.st);
clu_list = sp.cids(good_cells);
ACG_temp = zeros(numel(clu_list),lags+1,n_it);
PXX_temp = zeros(numel(clu_list),1025,n_it);
for iCell=1:nnz(good_cells)
    cid=clu_list(iCell);
    spike_id = sp.clu==cid;
    spike_t = sp.st(spike_id);
    quantiles = NaN(2,4001);
    
    max_t = spike_t(end);
    shuffles=max_t*rand(n_it,1);
    tmp_acg=zeros(n_it,lags+1);
    for num_it = 1:n_it
        shuffle_idx = shuffles(num_it);
        spike_t_shuffled = mod(spike_t+shuffle_idx,max_t)+abs_min;
        
        
        [~,~,spike_idx] = histcounts(spike_t_shuffled,post);
        spike_distance = distance(spike_idx);
        [aa,~]=histcounts(spike_distance,edges);
        
        %moothing
        firing_rate = aa./dT;
        firing_rate =zscore(firing_rate);
        [acg,spacing] = xcorr(firing_rate,lags,'coeff');
        ACG_temp(iCell,:,num_it)=acg((lags+1):end);
        pxx=pwelch(firing_rate,[],[],[],100/binsize);
        PXX_temp(iCell,:,num_it)=pxx;
    end
end

upper_bound = quantile(ACG_temp,0.99,3);
upper_bound_pxx = quantile(PXX_temp,0.99,3);
%%
cell_nr = 6;
figure
%errorbar(1:50,mean(ACG_temp(sid(cell_nr),1:50,:),3),quantile(ACG_temp(sid(cell_nr),1:50,:),.05,3)',quantile(ACG_temp(sid(cell_nr),1:50,:),.95,3))
hold on

plot(mean(ACG_temp(sid(cell_nr),1:50,:),3),'b')
plot(quantile(ACG_temp(sid(cell_nr),1:50,:),0.05,3),'b--')
plot(quantile(ACG_temp(sid(cell_nr),1:50,:),0.95,3),'b--')


plot(ACG(sid(cell_nr),1:50))
%%
figure
errorbar(1:50,mean(ACG_temp(sid(50),1:50,:),3),quantile(ACG_temp(sid(50),1:50,:),.05,3)',quantile(ACG_temp(sid(50),1:50,:),.95,3))