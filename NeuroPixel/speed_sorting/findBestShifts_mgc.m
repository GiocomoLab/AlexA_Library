function [data_out] = findBestShifts_mgc(data,ops)

%prepare variables
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
    depth = data.anatomy.tip_distance' - data.anatomy.z2;
end

reg = reg(data.sp.cgs==2);
if ~isempty(sub_reg)
    sub_reg=sub_reg(data.sp.cgs==2);
end
depth = depth(data.sp.cgs==2);
good_cells=data.sp.cids(data.sp.cgs==2);
factors = ops.factors;
trials = ops.trials;
nT = numel(trials);
edges = ops.edges;

    idx=triu(true(nT),1);


delays = ops.factors;
all_stability = nan(numel(good_cells),numel(delays));
st_orig = data.sp.st;
dat2=data;
 for dIdx = 1:numel(delays)
        
        dat2.sp.st = st_orig+delays(dIdx);
        frMat = calcTrialFRMat(good_cells,trials,dat2,ops);
        %corrMat = nan(numel(good_cells),numel(trials),numel(trials));
        for i = 1:numel(good_cells)
            fr_this = squeeze(frMat(i,:,:))';
            %fr_this = fr_this-nanmean(fr_this); % subtract mean on each trial
            %xcorr_this = xcorr(fr_this,0,'coeff')';
            %xcorr_this = reshape(xcorr_this,numel(trials),numel(trials));
            xcorr_this = corr(fr_this);
            all_stability(i,dIdx)=nanmean(xcorr_this(idx));
            %corrMat(i,:,:) = xcorr_this+repmat(diag(nan(numel(trials),1)),[1 1]); % make diagonals nan
        end
    end

good_idx = ismember(data.sp.clu,good_cells);
clu_tmp = data.sp.clu(good_idx);
st_tmp = data.sp.st(good_idx);
[a,~,clus]=unique(clu_tmp);
nClu = numel(a);



data_out.all_stability=all_stability;
data_out.region = reg;
data_out.sub_reg = sub_reg;
data_out.depth = depth;
data_out.CID = good_cells;
end

