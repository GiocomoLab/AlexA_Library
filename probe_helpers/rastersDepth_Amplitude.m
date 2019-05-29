function rastersDepth_Amplitude(vr_data,ks_metrics,image_save_dir)
params.TrackStart = 0;
params.SpatialBin = 2;
params.TrackEnd = 400;

metrics = ks_metrics;
sp = vr_data.sp;
posx=vr_data.posx;
post=vr_data.post;
trial = vr_data.trial;

cells_to_plot = sp.cids(sp.cgs==2); % all good cells



idx_time= double(metrics.spike_times)/30000 > sp.vr_session_offset & double(metrics.spike_times)/30000<=(sp.vr_session_offset+max(post));

spikeDepths = metrics.spike_depth;
spikeAmplitude = metrics.spike_amplitude;
% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    idx2=idx_time & metrics.spike_clusters==cells_to_plot(k);
    spike_depth(k) = median(spikeDepths(idx2));
end
spikeDepths = spikeDepths(idx_time);
spikeAmplitude = spikeAmplitude(idx_time);


%% make raster plots for all cells
h = figure('Position',[100 100 430 500]); hold on;

for k = 1:numel(cells_to_plot)
    
    fprintf('\tcell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
    
    % get spike times and index into post
    spike_id=sp.clu==cells_to_plot(k);
    spike_t = sp.st(spike_id);
    [~,~,spike_idx] = histcounts(spike_t,post);
    
    % make figure
    features={'spikeDepths','spikeAmplitude'};
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
        title(features{iF})
        xlabel(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
        xticks(''); yticks('');
        
    end
    
    
    saveas(h,fullfile(image_save_dir,sprintf('%d_%d.png',k,round(spike_depth(k)))),'png');
    
    clf;
end
close(h);
end
