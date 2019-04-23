function fig = spatial_raster(sp,clusterID,post,posx,trial,ax,dot_col)
if nargin < 6
    fig = figure('Position',[100 100 160 500]); hold on;
    ax=axes();
    dot_col='k';
end

if nargin == 6
    dot_col='k';
end

%% make raster plots
    %fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==clusterID);
    [~,~,spike_idx] = histcounts(spike_t,post);

    % make figure
    
 
    % plot spike raster
    axes(ax);
    plot(posx(spike_idx),trial(spike_idx),[dot_col '.']);
    %xlim([params.TrackStart params.TrackEnd]);
    ylim([0 max(trial)+1]);
    xticks(''); yticks('');

end