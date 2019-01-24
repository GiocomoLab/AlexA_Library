% plot average A period stability vs anatomical location of 
% unit based on histology and depth on probe
% MGC 1/21/2019

% make sure paths are correct
%restoredefaultpath
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\functions'));
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\spikes'));

% where to find data.
% these are the sessions with histology
data_dir = 'C:\Users\giocomolab\Dropbox\Work\neuropixels\data\';
session_name = {'npG2_1211_gain_1',...
    'npG2_1212_gaincontrast_1',...
    'npG5_1207_gain_1',...
    'npG5_1210_gaincontrast_2'};

% some params
params = readtable('C:\Users\giocomolab\Dropbox\Work\neuropixels\UniversalParams.xlsx');
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;
bl_block_size = 10; % size of baseline trial block in which to compute stability
  
% values to compute for each sessions
stability_all = {};
spike_depth_all = {};


% iterate over sessions
for session_num = 1:4
    load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));
  

    
    % good cells
    good_cells = sp.cids(sp.cgs>=1);
    
    % baseline trials, split into blocks (of size set above) to account for
    % slow drift (look into this more later...)
    baseline_trials = find(trial_gain==1 & trial_contrast==100);
    baseline_trials = reshape(baseline_trials,bl_block_size,numel(baseline_trials)/bl_block_size)';

    % save stability and cross corr values;
    stability = nan(numel(good_cells),size(baseline_trials,1));
    for k = 1:numel(good_cells)

        fprintf('session %d/%d, cell %d/%d\n',session_num,numel(session_name),k,numel(good_cells));

        spike_t = sp.st(sp.clu==good_cells(k));

        for i = 1:size(baseline_trials,1)
            posx_this = posx(ismember(trial,baseline_trials(i,:)));
            post_this = post(ismember(trial,baseline_trials(i,:)));
            spike_t_this = spike_t(spike_t>=min(post_this) & spike_t<max(post_this));
            trial_this = trial(ismember(trial,baseline_trials(i,:)));
            [~,~,idx_this] = histcounts(spike_t_this,post_this);
            stability(k,i) = calculateSpatialStability_trial(idx_this,posx_this,trial_this-min(trial_this)+1,params);            
        end
    end
    
    stability_all{session_num} = stability;
    spike_depth_all{session_num} = spike_depth;
    
end



% % get full range of stability
% stability_concat = [];
% for i = 1:numel(session_name)
%     stability_concat = [stability_concat; nanmean(stability_all{i},2)];
% end
% min_stability = min(stability_concat);
% max_stability = max(stability_concat);
% 
% % plot stability (coded by color and symbol size) vs anatomical location
% figure; hold on;
% for i = 1:numel(session_name)
%     
%     stability = stability_all{i};
%     spike_depth = spike_depth_all{i};
%     
%     % average stability over trial blocks for each unit
%     mean_stability = nanmean(stability,2);
%     
%     % scale it for plotting
%     plot_scale = (mean_stability-min_stability)/(max_stability-min_stability);
%     
%     % get 3D anatomical locations of units for this session
%     pos3D = origin(i,:)+spike_depth.*unit_vector(i,:);
%     
%     % make 3D plot (maybe there's a better way than a for loop...)
%     for j = 1:numel(spike_depth)
%         plot3(pos3D(j,1),pos3D(j,2),pos3D(j,3),'ko',...
%             'MarkerSize',0.1+10*plot_scale(j));
%     end
% end
% ylim([0 700]);
% zlim([-2500 1000]);
% axis equal
% xlabel('ML'); ylabel('AP'); zlabel('DV (from MEC border)');
% view(3);
% %%