function [h1] = plotExampleCellShift(data,shift_data,cluID,ops,subplot_ax_spikes,subplot_ax_speed)


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


cellIDX = find(data.sp.cids==cluID);
shift_idx = find(shift_data.CID == cluID);

all_factors = nanmean(shift_data.all_factors);


h1.fig=figure('Position',[440   236   319   562],'Renderer','Painters');

this_factor = all_factors(shift_idx);
used_factors(cellIDX)=this_factor;
% extract spike times for this cell
spike_id=data.sp.clu==cluID;
spike_t = data.sp.st(spike_id);
% convert to VR idx
[~,~,spike_idx] = histcounts(spike_t,data.post);
posx=mod(data.posx,max(ops.edges));
col = zeros(numel(spike_idx),3);
for iS = 1:numel(spike_idx)
    cT = data.trial(spike_idx(iS));
    col(iS,:)=get_color(data.trial_gain(cT),data.trial_contrast(cT));
end


factors = [0, this_factor];



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




[sp_blpre,sidx_pre]=sort(trial_speed_bl_pre,'descend');
[sp_blpo,sidx_post]=sort(trial_speed_bl_post,'descend');
[sp_gain,sidx_gain]=sort(trial_speed_gain,'descend');
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

h1.ax(1) = subplot(4,2,1);
scatter(data.posx(spike_idx),tmp(spike_idx),15,col,'.');
xlim(xl); xline(maxLoc);
ylim([0 nT])
%set(gca,'YDir','reverse')
%set(gca,'PlotBoxAspectRatio',[1 .5 1])
%set(gca,'OuterPosition',[0,0,1,.5])
h1.ax(2)=subplot(4,2,2);
scatter(sp_blpre,ops.bl_pre,30,'k','.')
hold on
scatter(sp_blpo,ops.bl_post,30,'k','.')
scatter(sp_gain,ops.gain_trials,30,'b','.')
%set(gca,'YDir','reverse')
%set(gca,'OuterPosition',[0,0,1,.5])

set(gca,'PlotBoxAspectRatio',[.5 1 1])


if ~isempty(subplot_ax_spikes)
    axes(subplot_ax_spikes)
    scatter(data.posx(spike_idx),tmp(spike_idx),15,col,'.');
    xlim(xl); xline(maxLoc);
    ylim([0 nT])
    %set(gca,'YDir','reverse')
    %set(gca,'PlotBoxAspectRatio',[1 .5 1])
    %set(gca,'OuterPosition',[0,0,1,.5])
    xlabel('Position [cm]')
    box off
    axes(subplot_ax_speed)
    scatter(sp_blpre,ops.bl_pre,30,'k','.')
    hold on
    scatter(sp_blpo,ops.bl_post,30,'k','.')
    scatter(sp_gain,ops.gain_trials,30,get_color(0.8,100),'.')
    %set(gca,'YDir','reverse')
    %set(gca,'OuterPosition',[0,0,1,.5])
    xlabel('Speed [cm/s]')
    set(gca,'PlotBoxAspectRatio',[.5 1 1])
    box off
end
figure(h1.fig);


subplot(4,2,3)
%scatter(data.posx(spike_idx),trial_sorted(spike_idx),15,col(spike_idx,:),'.');
scatter(posxhat(spike_idx),tmp(spike_idx),15,col,'.');
xlim(xl); xline(maxLoc);
ylim([0 nT])
xlabel(sprintf('%.3f',this_factor));
%set(gca,'PlotBoxAspectRatio',[1 .5 1])

%set(gca,'YDir','reverse')
FRMAT_SORTED=FRMAT(tmp_sid,:,:,cellIDX);
subplot(4,2,5)
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
%colorbar('southoutside')
colorbar('east')
subplot(4,2,7)
imagesc(ops.midpoints,trials,squeeze(FRMAT_SORTED(:,:,2)));
set(gca,'Clim',cl);
xlim(xl); xline(maxLoc);
%colorbar('southoutside')
colorbar('east')
for tr = trials
    patch([-15 0 0 -15]+xl(1)+5,[tr tr tr+1 tr+1]-0.5,get_color(data.trial_gain(tr),100),...
        'EdgeColor',get_color(data.trial_gain(tr),100));
end
colormap(cbrewer('seq','BuPu',20))


end

