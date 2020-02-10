function [data_out] = findPeakAlignement(data,ops)

%prepare variables
if isfield(data.anatomy,'parent_shifted')
    reg = data.anatomy.parent_shifted;
    sub_reg = data.anatomy.region_shifted;
else
    reg = data.anatomy.cluster_parent;
    sub_reg = data.anatomy.cluster_region;
end
if iscolumn(reg)
    reg = reg';
    sub_reg = sub_reg';
end

good_cells=data.sp.cids(data.sp.cgs==2);
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


good_idx = ismember(data.sp.clu,good_cells);
clu_tmp = data.sp.clu(good_idx);
st_tmp = data.sp.st(good_idx);
[uClu,~,clus]=unique(clu_tmp);
nClu = numel(uClu);
reg = reg(ismember(data.sp.cids,uClu));
sub_reg = sub_reg(ismember(data.sp.cids,uClu));


[spMapBL]=shiftAllMapsByFactor(ops,clus,st_tmp,nClu,data.posx,data.post,trial_sorted,speed,0);


spMapGain = spMapBL(end-3:end,:,:);
spMapBL = spMapBL(1:(numel(ops.trials)-4),:,:);

%figure
include_idx = -20/ops.BinWidth:1:20/ops.BinWidth;

allSlow = nan(size(spMapBL,3),numel(include_idx));
allFast = allSlow;
allGain = allSlow;
stability = nan(size(spMapBL,3),1);
idx = triu(numel(ops.trials));
allSpikes = cell(1,nClu);
maxInd=nan(1,nClu);
post_valid = data.post( ismember(data.trial,[ops.trials]));
posx_valid = data.posx( ismember(data.trial,[ops.trials]));
trial_valid = trial_sorted( ismember(data.trial,[ops.trials]));
for iC = 1:size(spMapBL,3)
    %%% neends fixing for small bins
    cc=corr(spMapBL(:,:,iC)');
    stability(iC)=nanmean(cc(idx));
    avFR = mean(spMapBL(:,:,iC));
    [ma,mi]=max(avFR(ops.search_range));
    [ma2,mi2]=max(avFR([ min(ops.search_range)-1 max(ops.search_range)+1] ));
    mi= mi+min(ops.search_range)-1;
    if ma>ma2
        maxInd(iC)=ops.midpoints(mi);
        trial_speed = getSpeedAroundPoint(speed,data.posx,data.trial,ops,ops.midpoints(mi),ops.speedWindow);
        trial_speed=trial_speed(1:(numel(ops.trials)-4)); %exclude gain trials
        [~,sidx]=sort(trial_speed);
        sl=mean(spMapBL(sidx(1:2),:,iC))/ma;
        fa = mean(spMapBL(sidx(end-1:end),:,iC))/ma;
        ga = mean(spMapGain(:,:,iC))/ma';
        allSlow(iC,:)=sl(mi+include_idx );
        allFast(iC,:)=fa(mi+include_idx );
        allGain(iC,:)=ga(mi+include_idx );
        this_cell_st = data.sp.st(data.sp.clu==uClu(iC));
        this_cell_spikeIdx = discretize(this_cell_st,post_valid);
        
        nonzero_idx = ~isnan(this_cell_spikeIdx);
        spike_pos = posx_valid(this_cell_spikeIdx(nonzero_idx));
        spike_trial = trial_valid(this_cell_spikeIdx(nonzero_idx));
        tmp = [spike_pos, spike_trial];
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
data_out.allSpikes = allSpikes;
data_out.speed = speed;
data_out.max_ind = maxInd;
end



