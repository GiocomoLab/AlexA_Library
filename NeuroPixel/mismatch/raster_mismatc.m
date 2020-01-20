% plots rasters for all cells in all sessions
% colors spikes by gain and contrast
% MGC 2/28/2019

%% params
% make sure paths are correct

% some params
params = readtable('UniversalParams.xlsx');
save_figs = true;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
data_dir = 'F:\NP_DATA\CL_OL_Towers';
%session_name = {'AA5_190809_gain_1'};
session_name = {};
sn = dir(fullfile(data_dir,'*mismatch*.mat'));
for iS = 1:numel(sn)
    %if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback'))
        session_name{end+1}=sn(iS).name(1:end-4);
    %end
end
% all the values we ever use
gains_all = [0.8 0.7 0.6 0.5 0.2];
contrasts_all = [100 50 20 10 5 2 0];

%% iterate over sessions
for session_num = 1:numel(session_name)    
    % load data
    fprintf('session %d/%d: %s\n',session_num,numel(session_name),session_name{session_num});
    load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));
    
    cells_to_plot = sp.cids(sp.cgs==2); % all good cells
    
%     if isfield(anatomy,'parent_shifted')
%         region = anatomy.parent_shifted;
%     else
%         region = anatomy.cluster_parent;
%     end
    %reg = startsWith(region,'VISp');
    cells_to_plot = sp.cids(sp.cgs==2);
    
    % make image dir if it doesn't exist
    image_save_dir = strcat('F:\images\',...
        session_name{session_num},'\pretty_rastersMismatch\');
    if exist(image_save_dir,'dir')~=7
        mkdir(image_save_dir);
    end

    % compute some useful information (like spike depths)
    [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
        templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

   
    
    
    
    % get spike depths
    spike_depth = nan(numel(cells_to_plot),1);
    for k = 1:numel(cells_to_plot)
        spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
    end

    % sort cells_to_plot by spike_depth (descending)
    [spike_depth,sort_idx] = sort(spike_depth,'descend');
    cells_to_plot = cells_to_plot(sort_idx);
    
    % get the mismatch trials
    speed_t=0.05;
% figure('Name',filenames{iF});; plot(speed)
%
if size(mismatch_trigger,1) ~=1
    mismatch_trigger=mismatch_trigger';
end

if nnz(mismatch_trigger>.5)>nnz(mismatch_trigger<.5)
    warning('flipping mismatch trigger');
    all_mm_trigs=strfind(mismatch_trigger<0.1,[0 0 1 1])+2;
else
   all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
end
speed=true_speed';
run_periods=smooth(speed,25)>speed_t;
run_window=-30:30;
possibles=strfind(run_periods',ones(1,length(run_window)))+floor(.5*length(run_window));
mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
mm_trials = trial(mm_trigs);    
    %% make raster plots for all cells
    h = figure('Position',[100 100 160 500]); hold on;
    
    for k =1:numel(cells_to_plot)
        
        cla;
        fprintf('\tcell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

        % get spike times and index into post
        spike_t = sp.st(sp.clu==cells_to_plot(k));
        [~,~,spike_idx] = histcounts(spike_t,post);

        % make figure

        % plot spike raster
        % baseline trials for all contrasts
   
            plot(posx(spike_idx),trial(spike_idx),'k.');
      
        % gain trials
        for jj = 1:numel(mm_trials)
            keep = trial(spike_idx)==mm_trials(jj);
            plot(posx(spike_idx(keep)),trial(spike_idx(keep)),'r.');
        end
        xlim([params.TrackStart params.TrackEnd]);
        ylim([0 max(trial)+1]);
        title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
        xticks(''); yticks('');
        
        % save fig
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
            %saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');
        end
        cla;
%         catch ME
%             warning(sprintf('did not work for cell %d',k))
    
    %end
    
        
    end
end