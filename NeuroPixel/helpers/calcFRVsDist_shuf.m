function [fr,distbincent,total_dist] = calcFRVsDist_shuf(cell_id,trials,dat,opt)
% cell_id:  which cell ids to analyze
% trials:   which trials to analyze; ** should be contiguous **
% dat:      the data structure, e.g. dat = load('npI1_0418_gain_3.mat')
% opt:      options
% 
% shuffles each cells spike times by a different amount
% MGC 1/28/2019

if ~isfield(opt,'dark')
    opt.dark = false;
end

% extract data for given trials
posx = dat.posx(ismember(dat.trial,trials));
post = dat.post(ismember(dat.trial,trials));
trial = dat.trial(ismember(dat.trial,trials));

% total distance run
total_dist = posx + opt.track_length * (trial-min(trial));
total_dist = total_dist - total_dist(1); % start from zero
if isfield(opt,'distbinedges')
    distbinedges = opt.distbinedges;    
elseif opt.dark
    distbinedges = 0:opt.SpatialBin:max(total_dist);
else
    distbinedges = 0:opt.SpatialBin:numel(trials)*opt.track_length;
end
distbincent = distbinedges(1:end-1)+opt.SpatialBin/2;

% time per distance bin
timeperbin = histcounts(total_dist,distbinedges);
timeperbin = timeperbin * opt.TimeBin;

% firing rate matrix
fr = nan(numel(cell_id),numel(distbincent));

for i = 1:numel(cell_id)
    % get spike times for this cell
    spike_t = dat.sp.st(dat.sp.clu==cell_id(i));
    spike_t = spike_t(spike_t>=min(post) & spike_t<=max(post));
    
    % shuffle spike times
    spike_t_shuf = spike_t-min(post);
    spike_t_shuf = mod(spike_t_shuf+unifrnd(20,max(post)-min(post),1),max(post)-min(post));
    spike_t_shuf = spike_t_shuf+min(post);
    
    % compute distance-binned firing rate
    [~,~,spike_idx] = histcounts(spike_t_shuf,post);
    fr_this = histcounts(total_dist(spike_idx),distbinedges);
    fr_this = fr_this./timeperbin;

    % interpolate missing values
    if sum(isnan(fr_this))>0
        fr_this = interp1(find(~isnan(fr_this)),fr_this(~isnan(fr_this)),1:numel(fr_this));
    end

    % smooth firing rate
    fr_this = gauss_smoothing(fr_this,opt.smoothSigma_dist/opt.SpatialBin);
    
    fr(i,:) = fr_this;
end

end