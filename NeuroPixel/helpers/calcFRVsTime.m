function [fr,spikeCount] = calcFRVsTime(cell_id,dat,opt,trials)
% cell_id:  which cell ids to analyze
% trials:   which trials to analyze; ** should be contiguous **
% dat:      the data structure, e.g. dat = load('npI1_0418_gain_3.mat')
% opt:      options
% MGC 12/15/2019

if exist('trials','var')==1
    % extract data for given trials
    post = dat.post(ismember(dat.trial,trials));
else
    post = dat.post;
end

% firing rate matrix
fr = nan(numel(cell_id),numel(post)-1);
spikeCount = fr;
for i = 1:numel(cell_id)
    % get spike times for this cell
    spike_t = dat.sp.st(dat.sp.clu==cell_id(i));
    spike_t = spike_t(spike_t>=min(post) & spike_t<=max(post));   
    
    % compute distance-binned firing rate
    fr_this = histcounts(spike_t,post)/opt.TimeBin;
    
    % smooth firing rate
    fr_this = gauss_smoothing(fr_this,opt.smoothSigma_time/opt.TimeBin);
    
    if sum(isnan(fr_this))>0
        keyboard
    end
    spikeCount(i,:)=histcounts(spike_t,post);
    fr(i,:) = fr_this;
end

fr = [fr fr(:,end)]; % make right size
spikeCount = [spikeCount spikeCount(:,end)];
end