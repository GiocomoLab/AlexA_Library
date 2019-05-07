% plots rasters for all cells in all sessions
% colors spikes by gain and contrast
% MGC 2/28/2019

%% params
% make sure paths are correct
%restoredefaultpath
%addpath(genpath('C:\Users\malcolmc\Dropbox\Work\neuropixels\functions'));
%addpath(genpath('C:\Users\malcolmc\Dropbox\Work\neuropixels\spikes'));

% some params
params = readtable('UniversalParams.xlsx');
save_figs = true;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
if ispc()
    data_dir = 'Z:\giocomo\attialex\NP_DATA\';
else
    data_dir = fullfile('/oak/stanford/groups/giocomo','attialex','NP_DATA');
end
session_name = {'npI5_0417_baseline_1','npF4_1025_gaincontrast_2'};



%% iterate over sessions
for session_num = 1:numel(session_name)    
    % load data
    fprintf('session %d/%d: %s\n',session_num,numel(session_name),session_name{session_num});
    load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));
    cells_to_plot = sp.cids(sp.cgs==2); % all good cells
    
    % make image dir if it doesn't exist
    if ispc()
    image_save_dir = strcat('F:\images\',...
        session_name{session_num},'\rasters_amplitude\');
    else
        image_save_dir = fullfile('/oak/stanford/groups/giocomo','attialex','images',session_name{session_num},'rasters_amplitude');
    end
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
    
    %% make raster plots for all cells
    h = figure('Position',[100 100 190 500]); hold on;

    for k = 1:numel(cells_to_plot)

        fprintf('\tcell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

        % get spike times and index into post
        spike_id=sp.clu==cells_to_plot(k);
        spike_t = sp.st(spike_id);
        [~,~,spike_idx] = histcounts(spike_t,post);

        % make figure
features={'spikeAmps'};
for iF =1:length(features)
    eval(['cf = ' features{iF} ';'])
subplot(1,length(features),iF)
scatter(posx(spike_idx),trial(spike_idx),2,cf(spike_id))
colormap winter
if nnz(unique(cf(spike_id)))>10
set(gca,'CLim',prctile(cf(spike_id),[25 75]))
end

        
        xlim([params.TrackStart params.TrackEnd]);
        ylim([0 max(trial)+1]);
        title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
        xticks(''); yticks('');
        end
        
        % save fig
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
        end
        cla;
    end
end