function plotExampleCells(data,shift_data,region,sn,ops)


%%
trials = ops.trials;
nT=numel(trials);
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
if ~isfield(ops,'speed_filter')
    ops.speed_filter = ops.filter;
end
if ~isempty(ops.speed_filter)
    speed_raw = conv(speed_raw,ops.speed_filter,'same');
end




if isfield(data.anatomy,'parent_shifted')
    region_this = data.anatomy.parent_shifted;
else
    region_this = data.anatomy.cluster_parent;
end
if iscolumn(region_this)
    region_this = region_this';
end

good_idx = data.sp.cgs==2 & startsWith(region_this,region);
%%
tmp = shift_data.all_factors;
stab = shift_data.all_stability >.5;
tmp(~stab)=nan;
all_factors = nanmean(tmp);
valid = sum(~isnan(shift_data.all_factors),1)>2;
all_factors(~valid)=nan;
used_factors = nan(nnz(good_idx),1);
good_cells = data.sp.cids(good_idx);
FRMAT = zeros(nT,ops.nBins,2,numel(good_cells));
col = zeros(numel(data.posx),3);
for iS = 1:numel(data.trial)
    cT = data.trial(iS);
    col(iS,:)=get_color(data.trial_gain(cT),data.trial_contrast(cT));
end
factor_cutoff = prctile(all_factors,5);

fig=figure('Position',[440   236   319   562],'Renderer','Painters','Visible','off');

for cellIDX=1:numel(good_cells)
    shift_idx = find(shift_data.CID == good_cells(cellIDX));
    this_factor = all_factors(shift_idx);
    if ~(this_factor<factor_cutoff)
        continue
    end
    used_factors(cellIDX)=this_factor;
    % extract spike times for this cell
    spike_id=data.sp.clu==good_cells(cellIDX);
    spike_t = data.sp.st(spike_id);
    % convert to VR idx
    [~,~,spike_idx] = histcounts(spike_t,data.post);
    posx=mod(data.posx,max(ops.edges));
    
    
    
    factors = [0, this_factor];
    
    % for each shift, calculate a spatial firing rate map and calculate
    % trial by trial correlation
    % average of this correlation is the STAB for this factor
    
    for iFactor = 1:numel(factors)
        spMatHat = zeros(nT,ops.nBins);
        posxhat = posx+factors(iFactor)*speed_raw;
        posxhat = mod(posxhat,max(ops.edges));
        spike_loc_hat = discretize(posxhat,ops.edges);
        
        
        occ_this = zeros(nT,ops.nBins);
        for iT=1:numel(spike_loc_hat)
            r=trial_sorted(iT);
            c=spike_loc_hat(iT);
            if ~isnan(r)
                occ_this(r,c)=occ_this(r,c)+1;
            end
        end
        occ_this = occ_this*ops.TimeBin;
        
        
        for ii=1:numel(spike_idx)
            if ~isnan(trial_sorted(spike_idx(ii)))
                spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))=spMatHat(trial_sorted(spike_idx(ii)),spike_loc_hat(spike_idx(ii)))+1;
            end
        end
        %spMatHat = medfilt1(spMatHat);
        %divide by occupancy
        spMatHat = spMatHat./occ_this;
        spMatHat = fillmissing(spMatHat,'pchip',2);
        iidx = (size(spMatHat,2)+1):(2*size(spMatHat,2));
        
        if ~isempty(ops.filter)
            spF = [spMatHat spMatHat spMatHat];
            spF = convn(spF,ops.filter,'same');
            spMatHat = spF(:,iidx);
        end
        FRMAT(:,:,iFactor,cellIDX)=spMatHat;
        
    end
    
    tmp_c = corr(squeeze(FRMAT(ops.bl_pre,:,1,cellIDX))');
    tmp_idx = triu(true(numel(ops.bl_pre)),1);
    stab = nanmean(tmp_c(tmp_idx));
%     if ~(stab>.4)
%         continue
%     end
%     
    meanFR = mean(FRMAT(ops.bl_pre,:,1,cellIDX));
    [tmp_i]=discretize([60,360],ops.edges);
    sta=tmp_i(1);
    sto=tmp_i(2);
    [ma,mi]=max(meanFR(sta:sto));
    mi = mi+sta-1;
    maxLoc = ops.midpoints(mi);
    ops_temp = ops;
    ops_temp.trials = trials;
    trial_speed = getSpeedAroundPoint(speed_raw,data.posx,data.trial,ops_temp,maxLoc,ops_temp.speedWindow);
    trial_speed_bl_pre = trial_speed(ops.bl_pre);
    trial_speed_gain = trial_speed(ops.gain_trials);
    trial_speed_bl_post = trial_speed(ops.bl_post);
    
    
    
    
    [~,sidx_pre]=sort(trial_speed_bl_pre,'descend');
    [~,sidx_post]=sort(trial_speed_bl_post,'descend');
    [~,sidx_gain]=sort(trial_speed_gain,'descend');
    tmp_rank = 1:numel(ops.bl_pre);
    tmp_rank(sidx_pre)=tmp_rank;
    tmp_rank_gain = ops.gain_trials;
    tmp_rank_gain(sidx_gain)=tmp_rank_gain;
    tmp_rank_post = ops.bl_post;
    tmp_rank_post(sidx_post)=tmp_rank_post;
    tmp_rank = [tmp_rank, tmp_rank_gain, tmp_rank_post];
    tmp_sid = [sidx_pre, sidx_gain+max(ops.bl_pre), sidx_post+max(ops.gain_trials)];
    
    tmp=trial_sorted;
    for iT=1:numel(tmp)
        if ~isnan(tmp(iT))
            tmp(iT)=tmp_rank(tmp(iT));
        end
    end
    xl = [maxLoc-60 maxLoc+40];
    
    subplot(4,1,1)
    scatter(data.posx(spike_idx),tmp(spike_idx),15,col(spike_idx,:),'.');
    %xlim(xl); 
    xline(maxLoc);
    ylim([0 nT])
    set(gca,'YDir','reverse')
    
    subplot(4,1,2)
    scatter(data.posx(spike_idx),trial_sorted(spike_idx),15,col(spike_idx,:),'.');
    %scatter(posxhat(spike_idx),tmp(spike_idx),15,col(spike_idx,:),'.');
    %xlim(xl); 
    xline(maxLoc);
    ylim([0 nT])
    
    set(gca,'YDir','reverse')
    FRMAT_SORTED=FRMAT(tmp_sid,:,:,cellIDX);
    subplot(4,1,3)
    imagesc(ops.midpoints,trials,squeeze(FRMAT_SORTED(:,:,1)));
    cl = get(gca,'Clim');
    cl = floor(cl);
    cl(2)=cl(2)+1;
    set(gca,'Clim');
    xlim(xl); xline(maxLoc);
    
    for tr = trials
        patch([-15 0 0 -15]+xl(1)+5,[tr tr tr+1 tr+1]-0.5,get_color(data.trial_gain(tr),100),...
            'EdgeColor',get_color(data.trial_gain(tr),100));
    end
        colorbar('eastoutside')

    subplot(4,1,4)
    imagesc(ops.midpoints,trials,squeeze(FRMAT_SORTED(:,:,2)));
    set(gca,'Clim',cl);
    xlabel(sprintf('%.3f',used_factors(cellIDX)))
    xlim(xl); xline(maxLoc);
    colorbar('eastoutside')
    for tr = trials
        patch([-15 0 0 -15]+xl(1)+5,[tr tr tr+1 tr+1]-0.5,get_color(data.trial_gain(tr),100),...
            'EdgeColor',get_color(data.trial_gain(tr),100));
    end
    savename=sprintf('%s_c%d_t%d:%d.pdf',sn,good_cells(cellIDX),trials(1),trials(end));
    saveas(gcf,fullfile(ops.savepath,savename))
    %pause
    clf
end
close(fig)
end