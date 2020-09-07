
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
%%
spike_idx = sp.clu==391;

[a,b]=discretize(sp.st(spike_idx),post);
features={'spikeAmps','spikeDepths'};
for iF =1:length(features)
    eval(['cf = ' features{iF} ';'])
figure;
scatter(posx(a),trial(a),2,cf(spike_idx))
if nnz(unique(cf(spike_idx)))>10
set(gca,'CLim',prctile(cf(spike_idx),[25 75]))
end
end

%%
% amplitudes = zeros(length(good_cells),max(trial));
% for iT=1:max(trial)
%     start = min(post(trial==iT));
%     stop = max(post(trial==iT));
%     for iC=1:length(good_cells)
%     spike_idx = sp.st>=start & sp.st<=stop & sp.clu==good_cells(iC);
%     amp = median(spikeAmps(spike_idx));
%     amplitudes(iC,iT) = amp;
%     end
% end

%%
figure
for iC=1:length(cells_to_plot)

tmp=correlation_All(:,:,iC);
imagesc(tmp)
set(gca,'CLim',prctile(tmp(:),[20 80]))
yyaxis right
plot(amplitudes(iC,:),'k','LineWidth',2)
ylim([min(amplitudes(iC,:)) max(amplitudes(iC,:))+1*range(amplitudes(iC,:))]);
pause
clf
end
%%
