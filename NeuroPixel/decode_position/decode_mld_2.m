%% baseline decoder

ops = load_default_opt;
ops.xbin = 2;
ops.xbinedges = 0:ops.xbin:400;
ops.xbincent = ops.xbinedges(1:end-1)+ops.xbin/2;
ops.nBins = numel(ops.xbincent);
    lw=ceil(3*ops.smoothSigma_dist);
    wx=-lw:lw;
    gw=exp(-wx.^2/(2*ops.smoothSigma_dist^2)); gw=gw/sum(gw);
ops.edges = ops.xbinedges;
% ops.filter = gausswin(11);
ops.filter = gw';
ops.trials=5:20;
ops.num_pcs_decoding = 10;
regions = {'MEC','VISp','RS'};
%%
%root = '/oak/stanford/groups/giocomo';
root = '/Volumes/Samsung_T5';
data_path = fullfile(root,'/attialex/NP_DATA/');
savepath = fullfile(root,'/attialex/mld_classifier2/');
if ~isfolder(savepath)
    mkdir(savepath)
end
files = dir(fullfile(data_path,'*.mat'));
animal_list = {};
valid_files = {};
for iF=1:numel(files)
    [~,sn]=fileparts(files(iF).name);
    
    parts = strsplit(sn,'_');
    animal = parts{1};
    date = parts{2};
    animal_date = [animal, '_',date];
    if ~contains(sn,{'mismatch','playback','dark'}) && ~ismember(animal_date,animal_list)
        valid_files{end+1}=files(iF).name;
        animal_list{end+1}=animal_date;
    end
end
%%
% p=gcp('nocreate');
% if isempty(p)
%     parpool(12);
% end

%%
for iF=1:numel(valid_files)
    data = load(fullfile(data_path,valid_files{iF}));
    [~,sn]=fileparts(valid_files{iF});
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
    
    good_cells_idx = data.sp.cgs == 2 & startsWith(reg,regions);
    if nnz(good_cells_idx)<2
        continue
    end
    if ~all(ismember(ops.trials,data.trial))
        continue
    end
    
    tmp_reg = reg(good_cells_idx);
    good_cells = data.sp.cids(good_cells_idx);
    good_idx = ismember(data.sp.clu,good_cells);
    clu_tmp = data.sp.clu(good_idx);
    st_tmp = data.sp.st(good_idx);
    
    [a,~,clus]=unique(clu_tmp);
    if isfield(data.anatomy,'depth')
        depth = data.anatomy.depth(good_cells_idx);
    else
        if ~isnan(data.anatomy.z2)
            depth = data.anatomy.tip_distance(good_cells_idx)-data.anatomy.z2;
        else
            depth = data.anatomy.tip_distance(good_cells_idx)-1000;
        end
    end
    if ~isrow(depth)
        depth=depth';
    end
    
    nClu = numel(a);
    ops_temp = ops;
    ops_temp.trials = 1:max(data.trial);
    [spMap,firingRate] = shiftAllMapsByFactor(ops_temp,clus,st_tmp,nClu,data.posx,data.post,data.trial,0,0);
    maxlag=100;
    XC=zeros(size(spMap,3),2*maxlag+1);
    stab = nan(1,size(spMap,3));
    for cellIDX =1:size(spMap,3)
        t1=squeeze(spMap(5:20,:,cellIDX));
        cc = corr(t1');
        c_idx = triu(true(size(t1,1)),1);
        stab(cellIDX) = nanmean(cc(c_idx));
    end
    
    if nnz(stab>.5)<5
        continue
    end
    good_cells = good_cells(stab>.5);
    region = tmp_reg(stab>.5);
    spMap = spMap(ops.trials,:,stab>.5);
    % the next part is a bit cumbersome, first turning it into long matrix
    % to zscore, then back again
    spMapFlat = zeros(size(spMap,1)*size(spMap,2),size(spMap,3));
    for iT=1:size(spMap,1)
        idx = (1:200)+(iT-1)*200;
        spMapFlat(idx,:)=squeeze(spMap(iT,:,:));
    end
    spMapFlatZ = zscore(spMapFlat,0,1);
    spMapNew = zeros(size(spMap));
    for iT=1:size(spMap,1)
        idx = (1:200)+(iT-1)*200;
        spMapNew(iT,:,:)=spMapFlatZ(idx,:);
    end
   
    speed = calcSpeed(data.posx,ops);
    trial_idx = ismember(data.trial,ops.trials);
    idx = speed>ops.SpeedCutoff & trial_idx;
   
    
   
   
    [~,~,posbin] = histcounts(posx,ops.xbinedges);
    posbin(posbin==0) = 1;
    pred_pos = [];
    for iFold = 1:numel(ops.trials)
        take_idx = true(1,numel(ops.trials));
        take_idx(iFold)=false;
     
        tc = squeeze(mean(spMapNew(take_idx,:,:),1)); % all trials except the one to be tested
        
        dot_prod = tc * squeeze(spMapNew(iFold,:,:))';
        [~,max_bin] = max(dot_prod);
        pred_pos = cat(2,pred_pos,ops.xbincent(max_bin));
        
    end

    mf = matfile(fullfile(savepath,sn),'Writable',true);
    mf.cluID = good_cells;
    mf.region = region;
    mf.pred_pos = pred_pos;
    mf.spMap = spMapNew;
end