% plots rasters for all cells in a session
% in one giant plot
% MGC 2/26/2019
session_name = 'npH3_0404_mismatch_1';
image_save_dir = strcat('F:\Data\npAna\images\',session_name,'\pretty_rasters\');
% make image dir if it doesn't exist
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end

%% analysis
% load data
load(fullfile(data_dir,strcat(session_name,'.mat')));
%%
cells_to_plot = sp.cids(sp.cgs==2); % all good cells

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


for k = 1:numel(cells_to_plot)

    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
    fig = spatial_raster(sp,cells_to_plot(k),post,posx,trial);
    offset=max(trial);
    fig=spatial_raster(open_loop.sp,cells_to_plot(k),open_loop.post,open_loop.posx,open_loop.trial+offset,fig,'r');
    % get spike times and index into post
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    xlim([round(min(posx)) round(max(posx))])
    % save fig

    saveas(fig,fullfile(image_save_dir,sprintf('%d.png',k)),'png');

    close;
end