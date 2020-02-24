%%
ops.region = 'MEC';
ops.binsize=5;
ops.smoothSigma = 4;
ops.maxlag = 10;
ops.filter = 7;
ops.nIterations = 300;
ops.quantile = .99;
%f_vec = linspace(0,1/ops.binsize/2,500);
p=700:-1:10;
f_vec=[0 1./(800:-.5:10)];
%sigma =1.5;
%window=floor(sigma*5/2)*2+1;
%fi = fspecial('gaussian',[1 window],sigma);.
fi=gausswin(ops.filter)';
ops.fi=fi/sum(fi);
% [filenames,triggers] = getFilesCriteria(ops.region,0,0,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
savepath = '/oak/stanford/groups/giocomo/attialex/distance_codingfft2';
%savepath = 'F:/temp/xcorrFFT';
if ~isfolder(savepath)
    mkdir(fullfile(savepath,'data'))
    mkdir(fullfile(savepath,'images'))
end

%save(sprintf('%s/ops.mat',savepath),'ops')
%%
files = dir('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*dark*');
%files = dir('F:/NP_DATA/np*dark*');
filenames={};
for ii=1:numel(files)
    filenames{ii}=fullfile(files(ii).folder,files(ii).name);
end
%%
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
        [~,sn,~]=fileparts(filenames{iF});

        datafile = sprintf('%s/data/%s.mat',savepath,sn);
        if isfile(datafile)
            disp(datafile)
            continue
        end
        
        speed = diff(data.posx);
        
        speed(speed<-10) = NaN;
        speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
        speed = [0;speed];
        distance = cumsum(speed);
        total_dist = data.posx+400*(data.trial-min(data.trial));
        total_dist = total_dist - total_dist(1);
        distance= total_dist;
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
        spatialMap = convn(spatialMap,ops.fi,'valid');
        
        
        
        
        %good_cells = data.sp.cgs==2 & strcmp(data.anatomy.cluster_parent,'MEC')';
        good_cells = data.sp.cgs==2;
        depth = data.anatomy.tip_distance(good_cells);
        [~,sid]=sort(depth,'descend');
        zSpatialMap = zscore(spatialMap(good_cells,:),0,[2]);
        clus = data.sp.cids(good_cells);
        firing_rate_units = zeros(1,numel(clus));
        max_time = max(data.sp.st);
        for ii=1:numel(clus)
            clu_id = clus(ii);
            firing_rate_units(ii)=sum(data.sp.clu==clu_id)/max_time;
        end
        
        
        
        
        lags = 200;
        nUnits = size(zSpatialMap,1);
        ACG = zeros(nUnits,lags+1);
        
        for ii=1:nUnits
            tmp=xcorr(zSpatialMap(ii,:),lags,'coeff');
            ACG(ii,:)=tmp((lags+1):end);
        end
        %PXX=pwelch(zSpatialMap',[],[],f_vec,1/binsize);
        
        Fs = 1/ops.binsize; % Sampling frequency                          
        L = floor(size(zSpatialMap,2)/2)*2; % Length of signal (making sure it's even)
        f = Fs*(0:(L/2))/L; % Define the frequency domain f
        X = zSpatialMap(:,1:L)';
        Y=fft(X);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1,:);
        P1(2:end-1,:) = 2*P1(2:end-1,:);
        PXX=interp1(f,P1,f_vec);

        
        n_spikes = numel(data.sp.st);
        nUnits = size(zSpatialMap,1);
        n_it = ops.nIterations;
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
            
            max_t = spike_t(end);
            shuffles=max_t*rand(n_it,1);
            tmp_acg=zeros(n_it,lags+1);
            for num_it = 1:n_it
                shuffle_idx = shuffles(num_it);
                spike_t_shuffled = mod(spike_t+shuffle_idx,max_t-abs_min)+abs_min;
                
                
                [~,~,spike_idx] = histcounts(spike_t_shuffled,data.post);
                spike_distance = distance(spike_idx);
                [aa,~]=histcounts(spike_distance,edges);
                
                %moothing
                firing_rate = aa./dT; %help
                firing_rate(isnan(firing_rate))=0;
                firing_rate = convn(firing_rate,ops.fi,'valid');
                firing_rate =zscore(firing_rate);
                [acg,spacing] = xcorr(firing_rate,lags,'coeff');
                ACG_temp(iCell,:,num_it)=acg((lags+1):end);
                
                X = firing_rate(1:L);
                Y=fft(X);
                P2 = abs(Y/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                pxx=interp1(f,P1,f_vec);
                PXX_temp(iCell,:,num_it)=pxx;
            end
        end
        
        upper_bound = quantile(ACG_temp,ops.quantile,3);
        upper_bound_pxx = quantile(PXX_temp,ops.quantile,3);
        
        fig = figure('visible','off');
        subplot(1,3,1)
        [~,sid]=sort(depth,'descend');
        imagesc(ACG(sid,:),[0 0.5])
        [mi,miii]=min(abs(depth(sid)-data.anatomy.z2));
        hold on
        plot([1 100],[miii miii],'k--')
        %yline(miii,'k-');
        hold on
        keep = false(numel(sid),1);
        peak_list = struct();
        peak_list(nUnits).peak_val=[];
        peak_list(nUnits).peak_loc=[];
        peak_list(nUnits).quantile = [];
        for ii=1:nUnits
            
            %[p,ip]=findpeaks(ACG(sid(ii),1:50),'SortStr','descend');
            [p,ip]=findpeaks(ACG(sid(ii),:),'MinPeakHeight',0);
            
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
        subplot(1,3,2)
        tmpACG = ACG(keep,:);
        [~,tmpsid] = sort(depth(keep),'descend');
        imagesc(tmpACG(tmpsid,:),[0 0.4])
        
        [mi,miii]=min(abs(sort(depth(keep),'descend')-data.anatomy.z2));
        hold on
        plot([1 100],[miii miii],'k--')
        
        subplot(1,3,3)
        v_idx =1./f_vec > 25 &  1./f_vec<600;
        [map,mip]=max(PXX(v_idx,:));
        offset = strfind(v_idx,[0 1]);
        mip=mip+offset;
        i_vec=-1:1;
        %mip = mip+10;
        keep = false(size(map));
        for iC=1:numel(map)
            if all(PXX(mip(iC)+i_vec,iC)>upper_bound_pxx(iC,mip(iC)+i_vec)')
                keep(iC)=true;
            end
        end
        tmpACG = ACG(keep,:);
        [~,tmpsid] = sort(depth(keep),'descend');
        imagesc(tmpACG(tmpsid,:),[0 0.4])
        
        
        [~,sn,~]=fileparts(filenames{iF});
        saveas(fig,sprintf('%s/images/%s.png',savepath,sn))
        mec_depth=depth-data.anatomy.z2;
        m=matfile(sprintf('%s/data/%s',savepath,sn),'Writable',true);
        m.peak_list = peak_list;
        m.mec_depth = mec_depth;
        m.ACG=ACG;
        m.PXX=PXX;
	m.firing_rate = firing_rate_units;
        m.upper_bound_pxx=upper_bound_pxx;
        m.f_vec = f_vec;
        m.cluster_anatomy = data.anatomy.cluster_parent(good_cells);
        %save(sprintf('/oak/stanford/groups/giocomo/attialex/distance_coding/data/%s',sn),'peak_list','mec_depth','ACG','PXX','upper_bound_pxx','f_vec')
        
        close(fig)
    catch ME
        sprintf('in file %s \n',filenames{iF})
        [~,sn,~]=fileparts(filenames{iF});
        fi=fopen(sprintf('%s/data/%s.err',savepath,sn),'w');
        fprintf(fi,ME.message);
        fclose(fi);
        disp(ME.message)
    end
    
end
