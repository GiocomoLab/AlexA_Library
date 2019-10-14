% makes matrix of correlations between trials
% using single trial firing rates
% in range of gains data
%
% addition: 
% loads cell_info from file
% you can input which brain region you want to look at
% MGC 9/27/2019

%% make sure paths are correct

root_dir = '/home/users/attialex/';
%neuropix_folder = fullfile(root_dir,'Dropbox','Work','neuropixels');
addpath(genpath(fullfile(neuropix_folder,'AlexA_Library')));
addpath(genpath(fullfile(neuropix_folder,'spikes')));

data_dir = fullfile('/oak/stanford/groups/giocomo/attialex','NP_DATA');
image_save_dir = fullfile('/oak/stanford/groups/giocomo/attialex','images','peak_xcorr');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
cell_info_dir = fullfile(neuropix_folder,'matlab_scripts','cell_info');
cell_info_file = 'cell_info_ALL_Oct2019_v2';

%% params

% change these params
brain_region = 'VISp';
num_trials_total = 36; % plots trials 1 to num_trials_total
stab_thresh = 0.3; % minimum average stability within 4 trial blocks
maxlag = 100; % in cm, for finding peak lag across all cells

% these params mostly stay fixed
params = readtable('UniversalParams.xlsx'));
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;
track_length = params.TrackEnd-params.TrackStart;

% plot figs?
save_figs = true;
gains_plot = [1 0.8 0.7 0.6 0.5 0.2];
plot_colors_gain = [0 0 0; cool(4); 0 0 1];
plot_colors_trial_num = gray(num_trials_total);

%% load cell_info and session_name
idx = strcmp(cell_info.BrainRegion,brain_region);
session_name = unique(cell_info.Session(idx));

%% filter cells by brain region and stability
local_stab = cell_info.LocalStability;
numblocks = floor(num_trials_total/4);
local_stab = mean(local_stab(:,1:numblocks),2);
keep_cell = ismember(cell_info.Session,session_name) & ...
    strcmp(cell_info.BrainRegion,brain_region) & ...
    local_stab > stab_thresh;
cell_info = cell_info(keep_cell,:);
local_stab = local_stab(keep_cell);
session_name = unique(cell_info.Session);

%% iterate over sessions
peak_xcorr_all = nan(num_trials_total,num_trials_total,numel(session_name));
for session_num = 1:numel(session_name)
    % load data
    load(fullfile(data_dir,session_name{session_num}));
    trial_gain = trial_gain(1:num_trials_total);

    % total distance run since start of session
    total_dist = posx + track_length*(trial-1);
    
    % position bin edges
    posbinedges = 0:params.SpatialBin:track_length*num_trials_total;
    
    % calc running speed
    speed = calcSpeed(posx,params);
    
    % good cells
    good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{session_num}));
    
    % single trial firing rate matrix
    fr_mat = nan(numel(good_cells),numel(posbinedges)-1);
    for k = 1:numel(good_cells)
        fprintf('session %d/%d: %s, cell %d/%d\n',session_num,numel(session_name),...
            session_name{session_num},k,numel(good_cells));

        spike_t = sp.st(sp.clu==good_cells(k));
        [~,~,spike_idx] = histcounts(spike_t,post);
        
        % filter by speed (only count periods where animal is running)
        keep = speed > params.SpeedCutoff;
        post_filt = post(keep);
        total_dist_filt = total_dist(keep);
        spike_t_filt = spike_t(keep(spike_idx));
        spike_t_filt = spike_t_filt(spike_t_filt<=max(post_filt));
        [~,~,spike_idx_filt] = histcounts(spike_t_filt,post_filt);
        
        % now compute distance-binned firing rate
        timeperbin = histcounts(total_dist_filt,posbinedges);
        timeperbin = timeperbin * params.TimeBin;
        firing_rate = histcounts(total_dist_filt(spike_idx_filt),posbinedges);
        firing_rate = firing_rate./timeperbin;
        
        % interpolate missing values
        if sum(isnan(firing_rate))>0
            firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
        end

        % gaussian filter for smoothing
        smoothSigma = params.SmoothSigmaFR/params.SpatialBin;
        smoothWindow = floor(smoothSigma*5/2)*2+1;
        gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

        % smooth firing rate
        firing_rate_smooth = conv(firing_rate,gauss_filter,'same');
        
        % enter this cell into the fr matrix
        fr_mat(k,:) = firing_rate_smooth;
    end
    
    % reshape fr_mat into dims of cells x positions x trials
    fr_mat = reshape(fr_mat,numel(good_cells),numel(xbincent),num_trials_total);    
    
    % correlation matrix
    peak_xcorr_this = nan(size(fr_mat,3),size(fr_mat,3),maxlag/params.SpatialBin*2+1,...
        size(fr_mat,1));
    for i = 1:size(fr_mat,1)
        fprintf('computing xcorr mat: %d/%d\n',i,size(fr_mat,1));
        x = squeeze(fr_mat(i,:,:));
        x_demeaned = x-nanmean(x);
        y = xcorr(x_demeaned,maxlag/params.SpatialBin,'coeff')';
        y = reshape(y,num_trials_total,num_trials_total,maxlag/params.SpatialBin*2+1);
        peak_xcorr_this(:,:,:,i) = y;
    end
    
    % enter results into matrix of trials x trials x sessions
    % for each pair of trials, finds lag with maximum mean xcorr across cells, and uses that lag
    peak_xcorr_all(:,:,session_num) = nanmax(nanmean(peak_xcorr_this,4),[],3);
end
grand_mean = nanmean(peak_xcorr_all,3);

%% plot xcorr matrices (both avg and individual sessions)

% all sessions individually
image_save_dir_this = fullfile(image_save_dir,'individual_sessions',...
    strcat(session_file,'_',brain_region));
if exist(image_save_dir_this,'dir')~=7
    mkdir(image_save_dir_this);
end
for i = 1:size(peak_xcorr_all,3)
    x = squeeze(peak_xcorr_all(:,:,i));
    x = x-diag(diag(x));
    h=figure;
    imagesc(x);
    xlim([1 num_trials_total]);
    ylim([1 num_trials_total]);
    axis square;
    title(session_name{i},'Interpreter','none');
    colorbar;
    if save_figs
        saveas(h,fullfile(image_save_dir_this,session_name{i}),'png');
        saveas(h,fullfile(image_save_dir_this,session_name{i}),'pdf');
        saveas(h,fullfile(image_save_dir_this,session_name{i}),'fig');
    end
end

% average xcorr across sessions
x = grand_mean-diag(diag(grand_mean));
h=figure;
imagesc(x);
xlim([1 num_trials_total]);
ylim([1 num_trials_total]);
axis square;
title('avg over sessions');
colorbar;
if save_figs
    saveas(h,fullfile(image_save_dir,strcat(session_file,'_',brain_region)),'png');
    saveas(h,fullfile(image_save_dir,strcat(session_file,'_',brain_region)),'pdf');
    saveas(h,fullfile(image_save_dir,strcat(session_file,'_',brain_region)),'fig');
end

%% cluster trials

image_save_dir_this = fullfile(image_save_dir,'cluster_trials');
if exist(image_save_dir_this,'dir')~=7
    mkdir(image_save_dir_this);
end
h=figure;
ha=tight_subplot(3,1);
% cluster grand mean matrix and plot dendrogram
axes(ha(1));
dist_vec = squareform(1-grand_mean,'tovector');
tree = linkage(dist_vec,'average');
[~,~,outperm]=dendrogram(tree,0);
[~,trial_sort_idx] = sort(outperm);
xlim([0.5 num_trials_total+0.5]);
xticks(''); yticks('');
% plot trial info
axes(ha(2));
for i = 1:num_trials_total
    % what gain is this trial?
    gain_this = find(gains_plot==trial_gain(i));
    % sorted order
    tr = trial_sort_idx(i);
    % colored by gain val:
    patch([tr-0.5 tr-0.5 tr+0.5 tr+0.5],[0 0.5 0.5 0], plot_colors_gain(gain_this,:));
    % colored by trial num:
    patch([tr-0.5 tr-0.5 tr+0.5 tr+0.5],[0.5 1 1 0.5], plot_colors_trial_num(i,:));
end
xlim([0.5 num_trials_total+0.5]); ylim([0 1]);
xticks(''); yticks('');
% plot dist matrix
axes(ha(3));
x = grand_mean-diag(diag(grand_mean));
x_sorted = x(:,outperm);
x_sorted = x_sorted(outperm,:);
imagesc(x_sorted-diag(diag(x_sorted)));
xlim([0.5 num_trials_total+0.5]);
ylim([0.5 num_trials_total+0.5]);
xticks(''); yticks('');
if save_figs
    saveas(h,fullfile(image_save_dir_this,strcat(session_file,'_',brain_region)),'png');
    saveas(h,fullfile(image_save_dir_this,strcat(session_file,'_',brain_region)),'pdf');
    saveas(h,fullfile(image_save_dir_this,strcat(session_file,'_',brain_region)),'fig');
end