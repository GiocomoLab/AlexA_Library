%%load data
%dataset = load('F:\NP_DATA\npF4_1025_gaincontrast_2.mat');
function [aa_spikes,maxChan_spikes,spikefig,aa_spikes_MEC]=quantification_session(dataset,bins,acg_data)

if nargin==2
    drawCoherence = false;
else
    drawCoherence = true;
end

%% compute metrics (Depth)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(dataset.sp.temps, dataset.sp.winv, dataset.sp.ycoords, dataset.sp.spikeTemplates, dataset.sp.tempScalingAmps);
depth=zeros(length(dataset.sp.cgs),1);
for iC=1:length(dataset.sp.cgs)
    depth(iC)=mean(spikeDepths(dataset.sp.clu==dataset.sp.cids(iC)));
end
dataset.sp.spikeAmps = spikeAmps;
dataset.spikeDepths = spikeDepths;
%% generate a spike matrix that is n_depth bins by n_timebins
valid_idx = ismember(dataset.sp.clu,dataset.sp.cids(dataset.sp.cgs>1));

[~,mecTriggerChan] = min(abs(bins-(dataset.anatomy.z2-800))); % z2 is in distance from tip,take off 800 um 
[~,mecChan] = min(abs(bins-(dataset.anatomy.z2)));
discrete_depth = discretize(spikeDepths(valid_idx),bins); %depth bin for each spike bin 1 is at tip!!
spike_times = dataset.sp.st(valid_idx);
spike_times = double(spike_times);

time_bins = 0.002;
discrete_time = round(spike_times/time_bins)+1;
spikeMat = zeros(numel(bins),(ceil(max(spike_times))/time_bins)+1,'single');
% populate matrix
for iS=1:numel(spike_times)
    time_idx=discrete_time(iS);
    d_idx = discrete_depth(iS);
    
    if ~isnan(d_idx)
        spikeMat(d_idx,time_idx)=spikeMat(d_idx,time_idx)+1;
    end
    
end
%% calculate power spectrum for each depth

[PxxSpikes,FSpikes] = pwelch(spikeMat',[],[],0.5:0.5:100,1/time_bins);
% get power in theta range, normalized by other frequencies
theta_range=[4 12];
theta_idx = FSpikes>theta_range(1) & FSpikes<=theta_range(2);
rest_idx = ~theta_idx & FSpikes>1;
thetaPower = mean(PxxSpikes(theta_idx,:));
restPower = mean(PxxSpikes(rest_idx,:));
thetaPowerN = thetaPower./restPower;
%find depth with strongest theta power
[a,maxChan_spikes]=max(thetaPowerN);
%bandpass that vector around theta band
bp_spikes = bandpass(spikeMat(maxChan_spikes,:),theta_range,1/time_bins);
%find peaks
[pks,locs_spikes] = findpeaks(bp_spikes);
mec_spikes = bandpass(spikeMat(mecTriggerChan,:),theta_range,1/time_bins);
[~,locs_mec] = findpeaks(mec_spikes);

%% extract snippets around peaks in all channels
snpsSpikes=extract_snps(spikeMat,locs_spikes,'win',[-100 100]);
tvec_spikes = [-100:100]*time_bins;
aa_spikes=squeeze(mean(snpsSpikes,3));

snpsSpikesMEC = extract_snps(spikeMat,locs_mec,'win',[-100 100]);
aa_spikes_MEC = squeeze(mean(snpsSpikesMEC,3));
%% some basic plotting
aa_spikesNorm=bsxfun(@rdivide,aa_spikes,sum(aa_spikes,2));
aa_spikesNormMEC=bsxfun(@rdivide,aa_spikes_MEC,sum(aa_spikes_MEC,2));
%aa_spikes=aa_spikes-mean(aa_spikes,2);
spikefig=figure('Position',[113         558        1127         420],'visible','off');
subplot(1,4,[2 3])
%imagesc(flipud(spikeMat(:,1:2000)))
tmp = flipud(spikeMat(:,1:2000));
tmp(tmp==0)=nan;
h=imagesc(tmp);
set(h, 'AlphaData', ~isnan(tmp))
if ~isnan(dataset.anatomy.z2)
    %[~,mec_bin]=min(abs(bins-dataset.anatomy.z2+800));
    fstring = 'k--';
else
    [~,mecTriggerChan]=min(abs(bins-1000)); %because thats what I chose fo estimate MEC entry
    fstring = 'r--';   
end
hold on
yline(numel(bins)-mecTriggerChan,fstring);
yline(numel(bins)-maxChan_spikes,'g');
yline(numel(bins)-mecChan,'b')
legend({'800 b MEC','maxChan','MEC'})
set(gca,'YTick',linspace(0,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,2000,10),'XTickLabel',round(linspace(1*0.002,2000*0.002,11),2))

xlim([1000 2000])
%c = gray;
%c = flipud(c);
colormap(summer);

subplot(1,4,1)
imagesc(flipud(aa_spikesNorm))
set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
%yline(385-highest_channel,'k');
xline(numel(tvec_spikes)/2+.5)
title('triggered spikes')
[~,max_loc]=max(aa_spikes,[],2);
hold on
plot((max_loc(end:-1:1)),1:numel(max_loc),'ro')

if drawCoherence
    subplot(1,4,4)
    plot(acg_data.avg_corr,acg_data.corr_midpoints)
end




%%
