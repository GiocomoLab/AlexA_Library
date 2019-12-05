function plotRasterSpeedSort(data,params,savepath)

fig = figure('Position',[109   487   742   259],'visible','off'); hold on;
gains_all = [0.8 0.7 0.6 0.5 0.2];
contrasts_all = [100 50 20 10 5 2 0];
gains = sort(unique(data.trial_gain),'descend');
contrasts = sort(unique(data.trial_contrast),'descend');
gains = gains(2:end);

trial_sorted = data.trial_sorted;
region = data.region;
posWindow = data.posWindow;

[~,gain_plot_idx] = ismember(gains,gains_all);
plot_colors_gain = cool(numel(gains_all));

% get plot colors for contrasts

[~,contrast_plot_idx] = ismember(contrasts,contrasts_all);
plot_colors_contrast = gray(numel(contrasts_all)+1);
plot_colors_contrast = plot_colors_contrast(1:end-1,:);
good_cells = data.sp.cids(data.sp.cgs==2);

for cellIDX = 1:numel(good_cells)
    cluID = find(data.sp.cids==good_cells(cellIDX));

subplot(1,4,[1 2 3 4])
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
ylim([0 max(data.trial_sorted)+1]);
title(sprintf('c%d, %s',good_cells(cellIDX),region{cluID}));
xticks(''); yticks('');
for ij=[80 :80:320]
    xline(ij);
end
for ij=posWindow
    xline(ij,'r');
end


%
% subplot(1,4,4)
% plot(delay_list(bl_trials,2),delay_list(bl_trials,1),'.')
% axis square
% xlabel('speed diff')
% ylabel('delay')
% x1=min(delay_list(:,2));
% x2=max(delay_list(:,2));
% y1=[1 x1]*mod.Coefficients.Estimate;
% y2=[1 x2]*mod.Coefficients.Estimate;
% hold on
% plot([x1 x2],[y1 y2])

saveas(fig,fullfile(savepath,sprintf('%d.png',cellIDX)),'png');
clf
end
close(fig)
end