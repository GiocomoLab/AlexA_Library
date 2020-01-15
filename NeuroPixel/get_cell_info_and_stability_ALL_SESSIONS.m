% Makes a table of cell info for all sessions in the data directory.
% For each unit in each session:
%     - Assigns a unique ID
%     - Extracts anatomical brain region
%     - Computes spatial stability within blocks of 4 trials across entire
%       session ("Local Stability")
% MGC 9/18/2019     

% make sure paths are correct
restoredefaultpath
root_dir = '/home/users/attialex/';
%neuropix_folder = fullfile(root_dir,'Dropbox','Work','neuropixels');
addpath(genpath(fullfile(root_dir,'AlexA_Library')));
addpath(genpath(fullfile(root_dir,'spikes')));

% extract names of all sessions in data dir
data_dir = fullfile('/oak/stanford/groups/giocomo/attialex','NP_DATA');
session_name = dir(fullfile(data_dir,'np*_gain*.mat'));
session_name = {session_name.name}';
session_name = session_name(contains(session_name,'.mat'));
for i = 1:numel(session_name)
    session_name{i} = session_name{i}(1:end-4); % remove ".mat"
end

%% params

% change these params
blocksize = 4; % num trials in block for computing LocalStability

% these params mostly stay fixed
params = readtable(fullfile('UniversalParams.xlsx'));
gains_all = [1 0.8 0.7 0.6 0.5 0.2]; % all gains we ever use
contrasts_all = [100 50 20 10 5 2 0]; % all contrasts we ever use

%% create cell_info table 
cell_info = table();
cell_info.Session = repmat({''},1000*numel(session_name),1);
cell_info.CellID = nan(1000*numel(session_name),1);
cell_info.UniqueID = repmat({''},1000*numel(session_name),1);
cell_info.BrainRegion = repmat({''},1000*numel(session_name),1);
cell_info.LocalStability = nan(1000*numel(session_name),250);
cell_info.NumTrials = nan(1000*numel(session_name),1);
cell_info.Gains = repmat({0},1000*numel(session_name),1);
cell_info.Contrasts = repmat({0},1000*numel(session_name),1);
cell_info.SessionType = repmat({''},1000*numel(session_name),1);

%% iterate over sessions
counter_cell = 1;
for session_num = 1:numel(session_name)
    % load data
    clear anatomy;
    fprintf('\nloading data for session %d/%d: %s ...\n',...
        session_num,numel(session_name),session_name{session_num});
    load(fullfile(data_dir,session_name{session_num}));
    num_trials_total = max(trial);
    numblocks = floor(num_trials_total/blocksize);
    trialblock = ceil(trial/blocksize);
    
    % calc running speed
    speed = calcSpeed(posx,params);
    
    % good cells
    good_cells = sp.cids(sp.cgs==2)';
    
    % get anatomy info
    if exist('anatomy','var')
        if sum(contains(fieldnames(anatomy),'parent_shifted')) > 0
            brain_region = anatomy.parent_shifted(sp.cgs==2)';
        elseif sum(contains(fieldnames(anatomy),'cluster_parent')) > 0
            brain_region = anatomy.cluster_parent(sp.cgs==2)';
        else
            brain_region = {'missing'};
        end
    else
        brain_region = {'missing'};
    end
    
    % compute local stability
    stability_this = nan(numel(good_cells),250);
    for k = 1:numel(good_cells)
        fprintf('session %d/%d: %s, cell %d/%d: Computing local stability\n',...
            session_num,numel(session_name),session_name{session_num},k,numel(good_cells));
        spike_t = sp.st(sp.clu==good_cells(k));
        [~,~,spike_idx] = histcounts(spike_t,post);        
        for j = 1:numblocks
            % filter by current trial block and by speed
            % (only count periods where animal is running)
            keep = speed > params.SpeedCutoff & trialblock==j;
            posx_filt = posx(keep);
            post_filt = post(keep);
            trial_filt = trial(keep);
            spike_t_filt = spike_t(keep(spike_idx));
            spike_t_filt = spike_t_filt(spike_t_filt<=max(post_filt));
            [~,~,spike_idx_filt] = histcounts(spike_t_filt,post_filt);
            stability_this(k,j) = calculateSpatialStability_trial(spike_idx_filt,...
                posx_filt,trial_filt,params);
        end 
    end
    
    % enter data into big table of all cells from all sessions
    cell_info.Session(counter_cell:counter_cell+numel(good_cells)-1) = session_name(session_num);
    cell_info.CellID(counter_cell:counter_cell+numel(good_cells)-1) = good_cells;
    tmp = strsplit(session_name{session_num},'_');
    for k = 1:numel(good_cells)
        cell_info.UniqueID{counter_cell+k-1} = sprintf('%s_%s_c%d',tmp{1},tmp{2},good_cells(k));
    end
    cell_info.BrainRegion(counter_cell:counter_cell+numel(good_cells)-1) = brain_region;    
    cell_info.LocalStability(counter_cell:counter_cell+numel(good_cells)-1,:) = stability_this;
    cell_info.NumTrials(counter_cell:counter_cell+numel(good_cells)-1) = num_trials_total;
    % information on gain and contrast manips:
    gains = sort(unique(trial_gain),'descend');
    contrasts = sort(unique(trial_contrast),'descend');
    if sum(~ismember(gains,gains_all)) > 0
        gains = 1;
    end
    if sum(~ismember(contrasts,contrasts_all)) > 0
        contrasts = 100;
    end
    cell_info.Gains(counter_cell:counter_cell+numel(good_cells)-1) = repmat({gains},numel(good_cells),1);
    cell_info.Contrasts(counter_cell:counter_cell+numel(good_cells)-1) = repmat({contrasts},numel(good_cells),1);
    cell_info.SessionType(counter_cell:counter_cell+numel(good_cells)-1) = {''};
    counter_cell = counter_cell+numel(good_cells);
end
cell_info = cell_info(1:counter_cell-1,:);