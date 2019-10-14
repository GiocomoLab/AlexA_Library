% computes percent remapping for range of gains sessions
% modified 9/28/2019 MGC

%% make sure paths are correct
restoredefaultpath
tol = 0.0001; % for comparing floats
root_dir = '/home/users/attialex/';
%neuropix_folder = fullfile(root_dir,'Dropbox','Work','neuropixels');
addpath(genpath(fullfile(root_dir,'AlexA_Library')));
addpath(genpath(fullfile(root_dir,'spikes')));

data_dir = fullfile('/oak/stanford/groups/giocomo/attialex','NP_DATA');
image_save_dir = fullfile('/oak/stanford/groups/giocomo/attialex','images','remapping');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
cell_info_dir = fullfile('/oak/stanford/groups/giocomo/attialex');
cell_info_file = 'session_info_cell';
load(fullfile(cell_info_dir,cell_info_file));

%% params

% change these params
session_file = 'sessions_starting_with_small_gain_MEC_and_VISp';
brain_region = 'MEC';
num_trials_total = 34;
gains_to_analyze = 0.8;
stab_thresh = 0.3;
numreps_max = 1; % maximum number of reps for each gain change
max_lag = 100; % for finding xcorr peak of individual cells

% these params mostly stay fixed
params = readtable('UniversalParams.xlsx'));
numtrialspergainchange = 4;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;
track_length = params.TrackEnd-params.TrackStart;

% params for plotting
save_figs = true;

%% load cell_info and session_name
load(fullfile(cell_info_dir,cell_info_file));
load(fullfile(cell_info_dir,session_file));
local_stab = cell_info.LocalStability;
local_stab = nanmean(local_stab(:,1:floor(num_trials_total/numtrialspergainchange)),2);
keep_cell = ismember(cell_info.Session,session_name) & ...
    strcmp(cell_info.BrainRegion,brain_region) & ...
    local_stab > stab_thresh;
cell_info = cell_info(keep_cell,:);
local_stab = local_stab(keep_cell);
session_name = session_name(ismember(session_name,cell_info.Session));

%% iterate over sessions
% values to compute
% last entry = A1-A2 comparison or A-B comparison
crosscorr_all = nan(numel(session_name)*1000,numel(gains_to_analyze),numreps_max,numel(xbincent)*2-1,2);
counter_cell = 1;
for session_num = 1:numel(session_name)
    % load data
    load(fullfile(data_dir,session_name{session_num}));
    
    % calc running speed
    speed = calcSpeed(posx,params);

    % good cells
    good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{session_num}));

    % compute cross corr values
    crosscorr_this = nan(numel(good_cells),numel(gains_to_analyze),numreps_max,numel(xbincent)*2-1,2);
    for k = 1:numel(good_cells)
        fprintf('session %d/%d: %s, cell %d/%d\n',session_num,numel(session_name),...
            session_name{session_num},k,numel(good_cells));
        spike_t = sp.st(sp.clu==good_cells(k));
        for i = 1:numel(gains_to_analyze)
            gain_change_trials = find(trial_gain==gains_to_analyze(i));
            numreps_this = min(numreps_max,sum(abs(diff(trial_gain)+1-gains_to_analyze(i))<tol));
            counter_trial = 1;
            for j = 1:numreps_this
                trialblock = gain_change_trials(counter_trial)-2*numtrialspergainchange:...
                    gain_change_trials(counter_trial)+numtrialspergainchange-1;
                trialblock = reshape(trialblock,numtrialspergainchange,...
                    numel(trialblock)/numtrialspergainchange)';
                fr = nan(numel(trialblock)/numtrialspergainchange,numel(xbincent));
                
                for m = 1:numel(trialblock)/numtrialspergainchange
                    posx_this = posx(ismember(trial,trialblock(m,:)));
                    post_this = post(ismember(trial,trialblock(m,:)));
                    spike_t_this = spike_t(spike_t>=min(post_this) & spike_t<max(post_this));
                    trial_this = trial(ismember(trial,trialblock(m,:)));
                    [~,~,idx_this] = histcounts(spike_t_this,post_this);
                    fr(m,:) = calculateSmoothedFiringRate(idx_this,posx_this,params);
                end
                % A1-A2 xcorr
                crosscorr_this(k,i,j,:,1) = xcorr(fr(2,:)-mean(fr(2,:)),fr(1,:)-mean(fr(1,:)),'coeff');
                % A-B xcorr
                crosscorr_this(k,i,j,:,2) = xcorr(fr(3,:)-mean(fr(3,:)),fr(2,:)-mean(fr(2,:)),'coeff');
                counter_trial = counter_trial + numtrialspergainchange;
            end
        end       
    end
    
    crosscorr_all(counter_cell:counter_cell+numel(good_cells)-1,:,:,:,:) = crosscorr_this;
    counter_cell = counter_cell+numel(good_cells);
end
crosscorr_all = crosscorr_all(1:counter_cell-1,:,:,:,:);
% crosscorr_all = squeeze(crosscorr_all);

%% plots
xplot = -track_length+params.SpatialBin:params.SpatialBin:track_length-params.SpatialBin;
keep_xcorr = xplot>=-max_lag & xplot<=max_lag;
xplot = xplot(keep_xcorr);
plot_colors = [cool(4); 0 0 1];

%% avg xcorr curves (all sessions individually)
for i = 1:numel(session_name)
    figure; hold on;
    for j = 1:numel(gains_to_analyze)
        for k = 1:numreps_max
            keep_cell = strcmp(cell_info.Session,session_name{i});
            N_bl = sum(~isnan(crosscorr_all(keep_cell,j,k,1,1)));
            N_gc = sum(~isnan(crosscorr_all(keep_cell,j,k,1,2)));
            mean_bl = squeeze(nanmean(crosscorr_all(keep_cell,j,k,:,1)));
            mean_gc = squeeze(nanmean(crosscorr_all(keep_cell,j,k,:,2)));
            sem_bl = squeeze(nanstd(crosscorr_all(keep_cell,j,k,:,1)))/sqrt(N_bl);
            sem_gc = squeeze(nanstd(crosscorr_all(keep_cell,j,k,:,1)))/sqrt(N_gc);
            errorbar(xplot,mean_bl(keep_xcorr),sem_bl(keep_xcorr),'k-');
            errorbar(xplot,mean_gc(keep_xcorr),sem_gc(keep_xcorr),'-','Color',plot_colors(j,:));
            xlim([-100 100]);
            plot(xlim(),[0 0],'k--');
            plot([0 0],ylim(),'k--');
            title(sprintf('%s\nn = %d cells bl, %d cells gc',...
                strrep(session_name{i},'_','-'),N_bl,N_gc));
        end
    end
end

%% avg xcorr curve over all sessions
figure; hold on;
for j = 1:numel(gains_to_analyze)
    for k = 1:numreps_max
        N_bl = sum(~isnan(crosscorr_all(:,j,k,1,1)));
        N_gc = sum(~isnan(crosscorr_all(:,j,k,1,2)));
        mean_bl = squeeze(nanmean(crosscorr_all(:,j,k,:,1)));
        mean_gc = squeeze(nanmean(crosscorr_all(:,j,k,:,2)));
        sem_bl = squeeze(nanstd(crosscorr_all(:,j,k,:,1)))/sqrt(N_bl);
        sem_gc = squeeze(nanstd(crosscorr_all(:,j,k,:,1)))/sqrt(N_gc);
        errorbar(xplot,mean_bl(keep_xcorr),sem_bl(keep_xcorr),'k-');
        errorbar(xplot,mean_gc(keep_xcorr),sem_gc(keep_xcorr),'-','Color',plot_colors(j,:));
        xlim([-100 100]);
        plot(xlim(),[0 0],'k--');
        plot([0 0],ylim(),'k--');
        title(sprintf('%s (n=%d sessions)',brain_region,numel(session_name)));
        xlabel('cm')
        ylabel('xcorr');
    end
end

%% compute pct remapping + map shift from mean xcorr for each session
peak_bl = nan(numel(session_name),numel(gains_to_analyze),numreps_max,1);
peak_gc = nan(numel(session_name),numel(gains_to_analyze),numreps_max,1);
shift_bl = nan(numel(session_name),numel(gains_to_analyze),numreps_max,1);
shift_gc = nan(numel(session_name),numel(gains_to_analyze),numreps_max,1);
for i = 1:numel(session_name)
    keep_cell = strcmp(cell_info.Session,session_name{i});
    for j = 1:numel(gains_to_analyze)
        for k = 1:numreps_max
            mean_bl = squeeze(nanmean(crosscorr_all(keep_cell,j,k,:,1)));
            mean_gc = squeeze(nanmean(crosscorr_all(keep_cell,j,k,:,2)));
            [peak_bl(i,j,k),shift_bl(i,j,k)] = max(mean_bl);
            [peak_gc(i,j,k),shift_gc(i,j,k)] = max(mean_gc);
        end
    end
end
pct_remap = peak_gc./peak_bl;
shift_gc = (shift_gc-numel(xbincent))*params.SpatialBin;
shift_bl = (shift_bl-numel(xbincent))*params.SpatialBin;