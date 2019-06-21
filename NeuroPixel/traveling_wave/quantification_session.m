%% compute metrics (Depth)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
        templatePositionsAmplitudes(dataset.sp.temps, dataset.sp.winv, dataset.sp.ycoords, dataset.sp.spikeTemplates, dataset.sp.tempScalingAmps);
    depth=zeros(length(dataset.sp.cgs),1);
    for iC=1:length(dataset.sp.cgs)
        depth(iC)=mean(spikeDepths(dataset.sp.clu==dataset.sp.cids(iC)));
    end
%%
valid_idx = ismember(dataset.sp.clu,dataset.sp.cids(dataset.sp.cgs>0));
bins = 0:40:3840;

discrete_depth = discretize(spikeDepths(valid_idx),bins);
spike_times = dataset.sp.st(valid_idx);
spike_times = double(spike_times);

time_bins = 0.002; 
discrete_time = round(spike_times/time_bins)+1;
spikeMat = zeros(numel(bins),(ceil(max(spike_times))/time_bins)+1);
for iS=1:numel(spike_times)
    time_idx=discrete_time(iS);
    d_idx = discrete_depth(iS);
    
    if ~isnan(d_idx)
    spikeMat(d_idx,time_idx)=spikeMat(d_idx,time_idx)+1;
    end

end
%%

[PxxSpikes,FSpikes] = pwelch(spikeMat',[],[],[],1/time_bins);
    
    theta_range=[4 12];
    theta_idx = FSpikes>theta_range(1) & FSpikes<=theta_range(2);
    rest_idx = ~theta_idx;
thetaPower = mean(PxxSpikes(theta_idx,:));
restPower = mean(PxxSpikes(rest_idx,:));
thetaPowerN = thetaPower./restPower;

[a,maxChan_spikes]=max(thetaPowerN);
bp_spikes = bandpass(spikeMat(maxChan_spikes,:),theta_range,1/time_bins);
[pks,locs_spikes] = findpeaks(bp_spikes);

%%
snpsSpikes=extract_snps(spikeMat,locs_spikes,'win',[-100 100]);
tvec_spikes = [-100:100]*time_bins;
aa_spikes=squeeze(mean(snpsSpikes,3));
aa_spikesNorm=bsxfun(@rdivide,aa_spikes,sum(aa_spikes,2));
%aa_spikes=aa_spikes-mean(aa_spikes,2);
spikefig=figure('Position',[113         558        1127         420]);
subplot(1,4,[2 4])
imagesc(flipud(spikeMat(:,1:2000)))
hold on
plot(locs_spikes(1:100),size(spikeMat,1)-maxChan_spikes*ones(1,100),'ro')
set(gca,'YTick',linspace(0,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,2000,10),'XTickLabel',round(linspace(1*0.002,2000*0.002,11),2))

subplot(1,4,1)
imagesc(flipud(aa_spikesNorm),[.05 8])
set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
%yline(385-highest_channel,'k');
xline(numel(tvec_spikes)/2+.5)
title('triggered spikes')
[~,max_loc]=max(aa_spikes,[],2);
hold on
plot((max_loc(end:-1:1)),1:numel(max_loc),'ro')



%%
