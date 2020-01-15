%%
files = dir('F:/NP_DATA/np*baseline*');
binsize=10;
DEPTH = [];
MAXF=[];
PI=[];
for iF=1:numel(files)
    
    try
    data = load(fullfile(files(iF).folder,files(iF).name));
    if ~ismember('anatomy',fieldnames(data))
        continue
    end
    
    speed = diff(data.posx);

speed(speed<-10) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = [0;speed];
distance = cumsum(speed);
distance(distance<0)=0;
edges=[0:binsize:(max(distance)+binsize)];

distance_trial = ones(size(data.post));




    idxVR=true(size(distance_trial));
    t_time=data.post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=data.sp.st<stop & data.sp.st>=start;
    [spM, dT]=getSpikeMatPosition(data.sp.st(idxNP),data.sp.clu(idxNP),distance,data.post(idxVR),'edges',edges,'max_clust',max(data.sp.clu)+1);

spatialMap=bsxfun(@rdivide,spM(data.sp.cids+1,:,:),dT);
spatialMap(isnan(spatialMap))=0;





good_cells = data.sp.cgs==2 & strcmp(data.anatomy.cluster_parent,'MEC')';
depth = data.anatomy.tip_distance(good_cells);
[~,sid]=sort(depth,'descend');
zSpatialMap = zscore(spatialMap(good_cells,:),0,[2]);



DEPTH = cat(1,DEPTH,depth-data.anatomy.z2);

lags = 100;
nUnits = size(zSpatialMap,1);
ACG = zeros(nUnits,lags+1);

for ii=1:nUnits
    tmp=xcorr(zSpatialMap(ii,:),lags,'coeff');
    ACG(ii,:)=tmp((lags+1):end);
end


n_spikes = numel(data.sp.st);
nUnits = size(zSpatialMap,1);
lags=100;
n_it = 200;
abs_min = min(data.sp.st);
clu_list = data.sp.cids(good_cells);
ACG_temp = zeros(numel(clu_list),lags+1,n_it);
PXX_temp = zeros(numel(clu_list),1025,n_it);
for iCell=1:nnz(good_cells)
    cid=clu_list(iCell);
    spike_id = data.sp.clu==cid;
    spike_t = data.sp.st(spike_id);
    if numel(spike_t)<3
        continue
    end
    quantiles = NaN(2,4001);
    
    max_t = spike_t(end);
    shuffles=max_t*rand(n_it,1);
    tmp_acg=zeros(n_it,lags+1);
    for num_it = 1:n_it
        shuffle_idx = shuffles(num_it);
        spike_t_shuffled = mod(spike_t+shuffle_idx,max_t)+abs_min;
        
        
        [~,~,spike_idx] = histcounts(spike_t_shuffled,data.post);
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

fig = figure('visible','off');
subplot(1,2,1)
[~,sid]=sort(depth,'descend');
imagesc(ACG(sid,1:50),[0 0.5])
hold on
keep = false(numel(sid),1);
peak_list = struct();
peak_list(nUnits).peak_val=[];
peak_list(nUnits).peak_loc=[];
peak_list(nUnits).quantile = [];
for ii=1:nUnits

    %[p,ip]=findpeaks(ACG(sid(ii),1:50),'SortStr','descend');
    [p,ip]=findpeaks(ACG(sid(ii),1:50),'MinPeakHeight',0);
    
    if ~isempty(ip)
        peak_list(sid(ii)).peak_val = p;
        peak_list(sid(ii)).peak_loc = ip;
        peak_list(sid(ii)).quantile = upper_bound(sid(ii),ip);
    PI(end+1)=ip(1);
            if p(1)>upper_bound(sid(ii),ip(1))
    plot(ip(1),ii,'kx')
    keep(sid(ii))=true;
    else
        plot(ip(1),ii,'rx')
    end

    else
        PI(end+1)=nan;
    end
end
subplot(1,2,2)
tmpACG = ACG(keep,:);
[~,tmpsid] = sort(depth(keep),'descend');
imagesc(tmpACG(tmpsid,1:50),[0 0.4])
saveas(fig,sprintf('F:/temp/xcorr_image/%s.png',files(iF).name))
mec_depth=depth-data.anatomy.z2;
save(sprintf('F:/temp/peak_list/%s',files(iF).name),'peak_list','mec_depth','ACG')

close(fig)
    catch ME
        disp(ME.message)
    end

end