% plots rasters for all cells in a session
% in one giant plot
% MGC 2/26/2019

%% data and params
% make sure paths are correct
restoredefaultpath
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\functions'));
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\spikes'));

% some params
params = readtable('UniversalParams.xlsx');
save_figs = true;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
data_dir = 'C:\Users\giocomolab\Dropbox\Work\neuropixels\data\';
session_name = 'npH5_0326_gaincontrast_1';
image_save_dir = strcat('C:\Users\giocomolab\Dropbox\Work\neuropixels\images\',session_name,'\pretty_rasters\');
% make image dir if it doesn't exist
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end

%% analysis
% load data
load(fullfile(data_dir,strcat(session_name,'.mat')));
cells_to_plot = sp.cids(sp.cgs==2); % all good cells

% compute some useful information (like spike depths)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

% find gain onset trials
gains = sort(unique(trial_gain),'descend');
gains = gains(2:end);
trial_gain_onsets = trial_gain;
trial_gain_onsets(trial_gain==1) = nan;
for i = numel(trial_gain_onsets):-1:2
    if trial_gain_onsets(i) == trial_gain_onsets(i-1)
        trial_gain_onsets(i) = nan;
    end
end
num_onsets_per_gain = sum(trial_gain_onsets==gains(1));
num_trials_per_block = sum(trial_gain==gains(1))/num_onsets_per_gain;
gain_onset_trials = nan(num_onsets_per_gain,numel(gains));
for k = 1:numel(gains)
    gain_onset_trials(:,k) = find(trial_gain_onsets==gains(k)); % gain change onset trials
end 

% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
end

% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
cells_to_plot = cells_to_plot(sort_idx);

%% make raster plots
plot_colors = cool(numel(gains));
for k = 1:numel(cells_to_plot)

    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    [~,~,spike_idx] = histcounts(spike_t,post);

    % make figure
    h = figure('Position',[100 100 160 500]); hold on;
    
    % plot patches over gain change trials
    for kk = 1:numel(gains)
        for jj = 1:num_onsets_per_gain
            start_trial = gain_onset_trials(jj,kk);
            end_trial = start_trial + num_trials_per_block-1;
            patch([params.TrackStart params.TrackEnd params.TrackEnd params.TrackStart],...
                [start_trial start_trial end_trial end_trial],plot_colors(kk,:));
        end
    end
    
    % plot spike raster
    plot(posx(spike_idx),trial(spike_idx),'k.');
    xlim([params.TrackStart params.TrackEnd]);
    ylim([0 max(trial)+1]);
    title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
    xticks(''); yticks('');

    % save fig
    if save_figs
        saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
        saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');
    end
    close;
end