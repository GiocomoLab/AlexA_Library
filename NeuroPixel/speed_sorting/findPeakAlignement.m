function [data_out] = findPeakAlignement(data,ops)

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
else % for MEC cases
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



[speed,speed_raw]=calcSpeed(data.posx,ops);
if ~isempty(ops.filter)
    speed_raw = conv(speed_raw,ops.filter,'same');
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

ops_tmp = ops.ops_shifts;
ops_tmp.trials=ops.trials(1:end-4); %only use pre gain trials
[data_shift,fighandles] = findBestShifts(data,ops_tmp); %align posx
[~,mi]=max(data_shift.all_stability,[],2);
factors = ops.factors(mi);

[spMapBL]=shiftAllMapsByFactor(ops,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed,0);


spMapGain = spMapBL(end-3:end,:,:);
spMapBL = spMapBL(1:(numel(ops.trials)-4),:,:);

%figure
include_idx = -20/ops.BinWidth:1:20/ops.BinWidth;

allSlow = nan(size(spMapBL,3),numel(include_idx));
allFast = allSlow;
allGain = allSlow;
stability = nan(size(spMapBL,3),1);
%idx = triu(numel(ops.trials));
idx = (triu(true(numel(ops.trials)-4),1));
idx_sim = triu(true(10),1);
similarity = stability;
allSpikes = cell(1,nClu);
maxInd=nan(1,nClu);
post_valid = data.post( ismember(data.trial,[ops.trials]));
posx_valid = data.posx( ismember(data.trial,[ops.trials]));
trial_valid = trial_sorted( ismember(data.trial,[ops.trials]));
for iC = 1:size(spMapBL,3)
    %%% neends fixing for small bins
    cc=corr(spMapBL(:,:,iC)');
    tmp=sum(spMapBL(:,:,iC),2);
    if nnz(tmp==0)<=2 %more than 2 trials without spikes
        stability(iC)=nanmean(cc(idx));
    end
    
    spMap = squeeze(cat(1,spMapBL(end-5:end,:,iC),spMapGain(:,:,iC)));
    mS=squeeze(spMap)';
    mS=mS-mean(mS);
    [xcorr_this,~]=xcorr(mS,'coeff',10);
    [xcorr_this,max_idx] = max(xcorr_this,[],1); % take max over lag
    xcorr_this = reshape(xcorr_this,10,10);
    similarity(iC) = nanmean(xcorr_this(idx_sim));
    
    avFR = mean(spMapBL(:,:,iC));
    [ma,mi]=max(avFR(ops.search_range));
    [ma2,mi2]=max(avFR([ min(ops.search_range)-1 max(ops.search_range)+1] ));
    mi= mi+min(ops.search_range)-1;
    if ma>ma2
        maxInd(iC)=ops.midpoints(mi);
        trial_speed = getSpeedAroundPoint(speed,data.posx,data.trial,ops,ops.midpoints(mi),ops.speedWindow);
        trial_speed=trial_speed(1:(numel(ops.trials)-4)); %exclude gain trials
        if nnz(isnan(trial_speed))>0
            keyboard
            error('nan speed')
        end
        [~,sidx]=sort(trial_speed);
        
        sl=mean(spMapBL(sidx(1:2),:,iC))/ma;
        fa = mean(spMapBL(sidx(end-1:end),:,iC))/ma;
        ga = mean(spMapGain(:,:,iC))/ma';
        allSlow(iC,:)=sl(mi+include_idx );
        allFast(iC,:)=fa(mi+include_idx );
        allGain(iC,:)=ga(mi+include_idx );
        this_cell_st = data.sp.st(data.sp.clu==uClu(iC));
        this_cell_spikeIdx = discretize(this_cell_st,post_valid);
        posx_hat = data.posx+factors(iC)*speed;
        posx_hat_valid = posx_hat( ismember(data.trial,[ops.trials]));
        
        
        nonzero_idx = ~isnan(this_cell_spikeIdx);
        spike_pos = posx_valid(this_cell_spikeIdx(nonzero_idx));
        spike_pos_hat = posx_hat_valid(this_cell_spikeIdx(nonzero_idx));
        spike_trial = trial_valid(this_cell_spikeIdx(nonzero_idx));
        tmp = [spike_pos, spike_trial,spike_pos_hat];
        allSpikes{iC}=tmp;
        %         subplot(2,1,1)
        %         plot(sl/ma)
        %         hold on
        %         plot(ga/ma)
        %         plot(fa/ma);
        %         plot(avFR/ma)
        %         legend({'slow','fast','gain','all'})
        %         xline(mi);
        %         xlim([mi-20,mi+20])
        %         subplot(2,1,2)
        %         plot(spike_pos,spike_trial,'.')
        %         gain_idx = spike_trial>(numel(ops.trials)-4);
        %         hold on
        %         plot(spike_pos(gain_idx),spike_trial(gain_idx),'r.')
        %         xline(ops.midpoints(mi));
        %         xlim(ops.midpoints(mi) +[-20 20])
        %         pause
        %         clf
    end
end
data_out.all_slow = allSlow;
data_out.all_fast = allFast;
data_out.all_gain = allGain;
data_out.region = reg;
data_out.subregion = sub_reg;
data_out.stability = stability;
data_out.similarity = similarity;
data_out.allSpikes = allSpikes;
data_out.speed = [speed, data.posx,data.trial];
data_out.max_ind = maxInd;
data_out.factors = factors;
data_out.CID = uClu;
data_out.depth = depth;
end



