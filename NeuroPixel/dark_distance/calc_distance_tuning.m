%% script looping over dark sessions for computing autocorellograms
ops = load_mismatch_opt;
ops.max_lag_autocorr = 800; %in cm
ops.num_shuffles = 300;
ops.SpatialBin = 2;
ops.dark = true;
%OAK = '/Volumes/Samsung_T5/attialex';
OAK = 'F:\Alex\';
matfiles = dir(fullfile(OAK,'matfiles_new','*MMdark*'));
outpath = fullfile(OAK,'distance_tuning_xcorr_only_3');
if ~isfolder(outpath)
    mkdir(outpath)
end

valid_files = true(1,numel(matfiles));
for iF = 1:numel(matfiles)
    if ismember('reward',matfiles(iF).name)
        valid_files(iF)=false;
    end
end

matfiles= matfiles(valid_files);

%%
p = gcp('nocreate');
if isempty(p)
    p = parpool(4);
end

%%
parfor iF=1:numel(matfiles)
    try
    
    [~,sn]=fileparts(matfiles(iF).name);
    
    if isfile(fullfile(outpath,matfiles(iF).name))
        sprintf('%s existsts, skipping \n',matfiles(iF).name)
        continue
    end
    data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
%     if ~isfield(data,'anatomy')
%         continue
%     end
%     if isfield(data.anatomy,'parent_shifted')
%         reg = data.anatomy.parent_shifted;
%     else
%         reg = data.anatomy.cluster_parent;
%     end
%     if iscolumn(reg)
%         reg=reg';
%     end
    good_idx = data.sp.cgs==2;
    good_cells = data.sp.cids(good_idx);
    %good_cells = data.sp.rf_cluster.cluster_id(startsWith(data.sp.rf_cluster.group,'good'))
    %reg = reg(good_idx);
    %depth = data.anatomy.tip_distance(good_idx);
    %mec_entry = data.anatomy.z2;
    
    %firing rate vs distance
    FRDist = calcFRVsDist(good_cells,[],data,ops);
    %FRDistN = FRDist - nanmean(FRDist,2);
    FRDistN = zscore(FRDist,0,2);
    firing_rate_units = zeros(1,numel(good_cells));
    max_time = max(data.sp.st);
    for ii=1:numel(good_cells)
        clu_id = good_cells(ii);
        firing_rate_units(ii)=nnz(data.sp.clu==clu_id)/max_time;
    end
    
    max_lag = ops.max_lag_autocorr/ops.SpatialBin;
    %lags = -ops.max_lag_autocorr:ops.SpatialBin:ops.max_lag_autocorr;
    nUnits = numel(good_cells);
    ACG = zeros(nUnits,max_lag+1);
    %compute acg for each unit
    peak_all = nan(numel(good_cells),1);
    peak_loc_all = nan(numel(good_cells),1);
    peak_prom_all = nan(numel(good_cells),1);
    %peak_shuf = nan(numel(good_cells),ops.num_shuffles);
    xplot = 0:ops.SpatialBin:ops.max_lag_autocorr;

    for ii=1:nUnits
        tmp=xcorr(FRDistN(ii,:),max_lag,'coeff');
        xc=tmp((max_lag+1):end);
        ACG(ii,:)=xc;
        [peaks,locs,~,prominence] = findpeaks(xc);
        if ~isempty(peaks)
            [peak_all(ii),idx] = max(peaks);
            peak_loc_all(ii) = xplot(locs(idx));
            peak_prom_all(ii) = prominence(idx);
        end
    end
    

    
    
    %shuffled data
    ACG_shuffled = zeros(nUnits,max_lag+1,ops.num_shuffles);
    peak_shuf = nan(numel(good_cells),ops.num_shuffles);

    for iS=1:ops.num_shuffles
        FRDistShuffle = calcFRVsDist_shuf(good_cells,[],data,ops);
        %FRDistShuffleN = FRDistShuffle - nanmean(FRDistShuffle,2);
        FRDistShuffleN = zscore(FRDistShuffle,0,2);
        for ii=1:nUnits
            tmp=xcorr(FRDistShuffleN(ii,:),max_lag,'coeff');
            xc=tmp((max_lag+1):end);
            ACG_shuffled(ii,:,iS)=xc;
            peaks = findpeaks(xc);
            if ~isempty(peaks)
                peak_shuf(ii,iS) = max(peaks);
            end
            
        end
    end
    
     pval = nan(numel(good_cells),1);
    for cIdx = 1:numel(good_cells)
        pval(cIdx) = sum(peak_shuf(cIdx,:)>=peak_all(cIdx))/ops.num_shuffles;
    end
    pval(isnan(peak_all) | sum(isnan(peak_shuf),2)>ops.num_shuffles/3) = nan;
    
    % only keep top percentile of shuffled data
    acg_upper_limits = prctile(ACG_shuffled,[91:2:99],3);
    %save
    chan_number =data.sp.waveform_metrics.peak_channel(ismember(data.sp.waveform_metrics.cluster_id,good_cells));
    mf = matfile(fullfile(outpath,sn),'Writable',true);
    mf.firing_rate = firing_rate_units;
    %mf.depth = depth;
    mf.chan_number = chan_number;
    %mf.region = reg;
    mf.ACG = ACG;
    mf.acg_upper_limits = acg_upper_limits;
    mf.pval = pval;
    mf.peak_all = peak_all;
    mf.peak_loc_all = peak_loc_all;
    mf.peak_prom_all = peak_prom_all;
    mf.peak_shuf = peak_shuf;
    %mf.mec_entry = mec_entry;
    catch ME
        fprintf('file: %d, msg: %s \n',iF,ME.message)
    end
    
end

%%
