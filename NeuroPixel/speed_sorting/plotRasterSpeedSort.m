function trial_rank = plotRasterSpeedSort(data,params,savepath,clusters2show,trials2show,ax)

if ~isempty(trials2show)
    height = ceil(numel(trials2show)/16)*85;
else
    height = ceil(max(data.trial)/16)*85;
end
if isempty(ax)
    fig = figure('Position',[109   387   742   height],'visible','on'); hold on;
    ax=axes();
else
    axes(ax);
    fig = gcf;
end
gains_all = [0.8 0.7 0.6 0.5 0.2];
contrasts_all = [100 50 20 10 5 2 0];
gains = sort(unique(data.trial_gain),'descend');
contrasts = sort(unique(data.trial_contrast),'descend');
gains = gains(2:end);
posWindow = data.posWindow;

if isfield(data,'trial_sorted')
    trial_sorted = data.trial_sorted;
else
 
    speed = calcSpeed(data.posx,params);
    
    trial_speed = zeros(1,max(data.trial));
    for ii = 1:max(data.trial)
        idx = data.posx>posWindow(1) & data.posx<posWindow(2) & data.trial == ii;
        tmp = mean(speed(idx));
        
        trial_speed(ii)=tmp;
        
        
    end
    if isempty(trials2show)
        [~,sid]=sort(trial_speed);
        trials2show = 1:max(data.trial);
    else
        ts=trial_speed;
        ts(~ismember(1:numel(ts),trials2show))=nan;
        [~,sid]=sort(ts);
    end
        
    
    trial_sorted = data.trial;
    trial_rank = 1:max(data.trial);
    trial_rank(sid)=trial_rank;
    for ii=1:max(data.trial)
                    idx = data.trial==ii;
        if ismember(ii,trials2show)
            trial_sorted(idx)=trial_rank(ii);
        else
            trial_sorted(idx)=nan;
        end
    end
    
end
if isfield(data.anatomy,'parent_shifted')
region = data.anatomy.parent_shifted;
else
    region = data.anatomy.cluster_parent;
end

[~,gain_plot_idx] = ismember(gains,gains_all);
plot_colors_gain = cool(numel(gains_all));

% get plot colors for contrasts

[~,contrast_plot_idx] = ismember(contrasts,contrasts_all);
plot_colors_contrast = gray(numel(contrasts_all)+1);
plot_colors_contrast = plot_colors_contrast(1:end-1,:);
if isempty(clusters2show)
good_cells = data.sp.cids(data.sp.cgs==2);
else
    good_cells = clusters2show;
end

for cellIDX = 1:numel(good_cells)
    cluID = find(data.sp.cids==good_cells(cellIDX));
    
    hold on
    % get spike times and index into post
    spike_t = data.sp.st(data.sp.clu==good_cells(cellIDX));
    [~,~,spike_idx] = histcounts(spike_t,data.post);
    
    
    for j = 1:numel(contrasts)
        keep = data.trial_contrast(data.trial(spike_idx))==contrasts(j) & ...
            data.trial_gain(data.trial(spike_idx))==1;
        plot(data.posx(spike_idx(keep)),trial_sorted(spike_idx(keep)),'.','Color',plot_colors_contrast(contrast_plot_idx(j),:));
    end
    % gain trials
    for j = 1:numel(gains)
        keep = data.trial_gain(data.trial(spike_idx))==gains(j);
        plot(data.posx(spike_idx(keep)),trial_sorted(spike_idx(keep)),'.','Color',plot_colors_gain(gain_plot_idx(j),:));
    end
    xlim([params.TrackStart params.TrackEnd]);
    ylim([0 max(trial_sorted)+1]);
    %title(sprintf('c%d, %s',good_cells(cellIDX),region{cluID}));
    ylabel(sprintf('c%d, %s',good_cells(cellIDX),region{cluID}));
    xticks(''); yticks('');
    for ij=[80 :80:320]
        xline(ij);
    end
    for ij=posWindow
        xline(ij,'r');
    end
    
    if ~isempty(savepath)
        saveas(fig,fullfile(savepath,sprintf('%d.png',cellIDX)),'png');
        clf
    end
    
end
if ~isempty(savepath)
    close(fig)
end
end