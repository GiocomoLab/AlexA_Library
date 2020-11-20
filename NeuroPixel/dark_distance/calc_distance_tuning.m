%% script looping over dark sessions for computing autocorellograms
ops = load_default_opt;
ops.max_lag_autocorr = 800; %in cm
ops.f_vec=[0 1./(800:-1:10)]; %firing rates to evaluate fft
ops.num_shuffles = 300;
%OAK = '/Volumes/Samsung_T5/attialex';
OAK = '/oak/stanford/groups/giocomo/attialex';
matfiles = dir(fullfile(OAK,'NP_DATA_corrected','np*dark*'));
outpath = fullfile(OAK,'distance_tuning_xcorr_only');
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
    p = parpool(12);
end

%%
parfor iF=1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    
    if isfile(fullfile(outpath,matfiles(iF).name))
        sprintf('%s existsts, skipping \n',matfiles(iF).name)
        continue
    end
    data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    if ~isfield(data,'anatomy')
        continue
    end
    if isfield(data.anatomy,'parent_shifted')
        reg = data.anatomy.parent_shifted;
    else
        reg = data.anatomy.cluster_parent;
    end
    if iscolumn(reg)
        reg=reg';
    end
    good_idx = data.sp.cgs==2;
    good_cells = data.sp.cids(good_idx);
    reg = reg(good_idx);
    depth = data.anatomy.tip_distance(good_idx);
    mec_entry = data.anatomy.z2;
    
    %firing rate vs distance
    FRDist = calcFRVsDist(good_cells,[1:max(data.trial)-1],data,ops);
    FRDistN = FRDist - nanmean(FRDist,2);
    
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
    for ii=1:nUnits
        tmp=xcorr(FRDistN(ii,:),max_lag,'coeff');
        ACG(ii,:)=tmp((max_lag+1):end);
    end
    
    %shuffled data
    ACG_shuffled = zeros(nUnits,max_lag+1,ops.num_shuffles);
    for iS=1:ops.num_shuffles
        FRDistShuffle = calcFRVsDist_shuf(good_cells,[1:max(data.trial)-1],data,ops);
        FRDistShuffleN = FRDistShuffle - nanmean(FRDistShuffle,2);
        
        for ii=1:nUnits
            tmp=xcorr(FRDistShuffleN(ii,:),max_lag,'coeff');
            ACG_shuffled(ii,:,iS)=tmp((max_lag+1):end);
        end
    end
    % only keep top percentile of shuffled data
    acg_upper_limits = prctile(ACG_shuffled,[91:2:99],3)
    %save
    mf = matfile(fullfile(outpath,sn),'Writable',true);
    mf.firing_rate = firing_rate_units;
    mf.depth = depth;
    mf.region = reg;
    mf.ACG = ACG;
    mf.acg_upper_limits = acg_upper_limits;
    mf.mec_entry = mec_entry;
    
end

%%
