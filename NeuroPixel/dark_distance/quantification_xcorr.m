%%
ops.region = 'MEC';
ops.binsize=10;
ops.smoothSigma = 4;
ops.maxlag = 10;
f_vec = linspace(0,1/ops.binsize/2,500);
% [filenames,triggers] = getFilesCriteria(ops.region,0,0,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
%%
files = dir('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*dark*')
filenames={};
for ii=1:numel(files)
    filenames{ii}=fullfile(files(ii).folder,files(ii).name);
end
p=gcp('nocreate');
if isempty(p)
    parpool(6);
end

%%
binsize=ops.binsize;

parfor iF=1:numel(filenames)
    
    try
        data = load(filenames{iF});
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
        
        
        
        
        lags = 100;
        nUnits = size(zSpatialMap,1);
        ACG = zeros(nUnits,lags+1);
        
        for ii=1:nUnits
            tmp=xcorr(zSpatialMap(ii,:),lags,'coeff');
            ACG(ii,:)=tmp((lags+1):end);
        end
        PXX=pwelch(zSpatialMap',[],[],f_vec,1/binsize);
        
        n_spikes = numel(data.sp.st);
        nUnits = size(zSpatialMap,1);
        lags=100;
        n_it = 200;
        abs_min = min(data.sp.st);
        clu_list = data.sp.cids(good_cells);
        ACG_temp = zeros(numel(clu_list),lags+1,n_it);
        PXX_temp = zeros(numel(clu_list),numel(f_vec),n_it);
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
                pxx=pwelch(firing_rate,[],[],f_vec,1/binsize);
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
                if p(1)>upper_bound(sid(ii),ip(1))
                    plot(ip(1),ii,'kx')
                    keep(sid(ii))=true;
                else
                    plot(ip(1),ii,'rx')
                end
                
            else
            end
        end
        subplot(1,2,2)
        tmpACG = ACG(keep,:);
        [~,tmpsid] = sort(depth(keep),'descend');
        imagesc(tmpACG(tmpsid,1:50),[0 0.4])
        [~,sn,~]=fileparts(filenames{iF});
        saveas(fig,sprintf('/oak/stanford/groups/giocomo/attialex/distance_coding/images/%s.png',sn))
        mec_depth=depth-data.anatomy.z2;
        m=matfile(sprintf('/oak/stanford/groups/giocomo/attialex/distance_coding/data/%s',sn));
        m.peak_list = peak_list;
        m.mec_depth = mec_depth;
        m.ACG=ACG;
        m.PXX=PXX;
        m.upper_bound_pxx=upper_bound_pxx;
        m.f_vec = f_vec;
        %save(sprintf('/oak/stanford/groups/giocomo/attialex/distance_coding/data/%s',sn),'peak_list','mec_depth','ACG','PXX','upper_bound_pxx','f_vec')
        
        close(fig)
    catch ME
        disp(ME.message)
    end
    
end