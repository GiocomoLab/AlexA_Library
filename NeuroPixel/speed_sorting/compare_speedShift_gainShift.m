function data_out = compare_speedShift_gainShift(data,ops)
if isfield(data.anatomy,'parent_shifted')
    reg = data.anatomy.parent_shifted;
else
    reg = data.anatomy.cluster_parent;
end

if isfield(data.anatomy,'parent_shifted')
    sub_reg = data.anatomy.region_shifted;
elseif isfield(data.anatomy,'cluster_region')
    
    sub_reg = data.anatomy.cluster_region;
else
    sub_reg = {};
end

if iscolumn(reg)
    reg = reg';
    sub_reg = sub_reg';
end

if isfield(data.anatomy,'depth_shifted')
    depth = data.anatomy.depth_shifted;
elseif isfield(data.anatomy,'depth')
    depth = data.anatomy.depth;
end
if isfield(data.anatomy,'z2') % for MEC cases
    depth = data.anatomy.tip_distance - data.anatomy.z2;
end

trials = ops.trials;


%create trial map that only contains numbers for trials to be included
trialMap = nan(1,numel(data.trial_gain));

cntr = 1;
for iT =1:numel(data.trial_gain)
    if ismember(iT,trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end




%recreate data.trial map that resets numbers of included trials and sets
%all else to nan
trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end




good_cells = data.sp.cids(data.sp.cgs==2);
good_idx = ismember(data.sp.clu,good_cells);
clu_tmp = data.sp.clu(good_idx);
st_tmp = data.sp.st(good_idx);
[uClu,~,clus]=unique(clu_tmp);
nClu = numel(uClu);
reg = reg(ismember(data.sp.cids,uClu));
depth = depth(ismember(data.sp.cids,uClu));
if ~isempty(sub_reg)
    sub_reg = sub_reg(ismember(data.sp.cids,uClu));
end
if nClu ~= numel(good_cells) % in case there are 'good cells' that don't have a spike in this data strct
    only_good = ~ismember(good_cells,uClu);
    only_good = good_cells(only_good);
    for iC = 1:numel(only_good)
        idx =  data.sp.cids==only_good(iC);
        data.sp.cgs(idx)=1; %set to MUA;
    end
end

[speed,~]=calcSpeed(data.posx,ops);
ops_tmp = ops; 
ops_tmp.trials=ops.trials(ops.bl_pre); % to exclude gain trials
[data_shift,~] = findBestShifts(data,ops); %align posx
[~,mi]=max(data_shift.all_stability,[],2);
factors = ops.factors(mi);


[spMapBL]=shiftAllMapsByFactor(ops,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed,0);
spMapShifted = shiftAllMapsByFactor(ops,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed,factors);

[a,b]=spMapXcorr(spMapBL,10,ops.BinWidth);
[ahat,bhat]=spMapXcorr(spMapShifted,10,ops.BinWidth);
nC=size(spMapBL,3);
similarity = nan(nC,1);
idx_sim = triu(true(10),1);
for iC = 1:nC
    
    
    mS = squeeze(spMapBL(ops.similarity_trials,:,iC))';
    mS=mS-mean(mS);
    [xcorr_this,~]=xcorr(mS,'coeff',10);
    [xcorr_this,~] = max(xcorr_this,[],1); % take max over lag
    xcorr_this = reshape(xcorr_this,numel(ops.similarity_trials),numel(ops.similarity_trials));
    similarity(iC) = nanmean(xcorr_this(idx_sim));
end
data_out.similarity = similarity;
data_out.region = reg;
data_out.corrMat = a;
data_out.shiftMat = b;
data_out.corrMatHat = ahat;
data_out.shiftMatHat = bhat;
data_out.factors = factors;
data_out.CID = uClu;
data_out.depth = depth;
end

