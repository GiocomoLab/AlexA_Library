myKsDir = 'F:\J4\npJ4_0511_baseline_g0\npJ4_0511_baseline_g0_imec0';
channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
ks = Neuropixel.KiloSortDataset(myKsDir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

nClips = 10;
clipDur = 200; % seconds

% load nClips one-sec samples
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
skip_seconds = 400;
sampStarts = round(linspace(lfpFs*skip_seconds, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(lfpFs*clipDur);



mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
n=1;
thisDat = double(mmf.Data.x(:, (1:nClipSamps)+sampStarts(n)));
thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));
%%
L = size(thisDat',1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(thisDat',[],[],[],lfpFs);
%%
theta_range=[4 12];
thetaPower = mean(Pxx(F>theta_range(1) & F<=theta_range(2),:));
alphaPower = mean(Pxx(F>0.1 & F<=theta_range(1),:));
omegaPower = mean(Pxx(F>theta_range(2) & F<=20,:));
figure
hold on
plot(10*log10(thetaPower))
hold on
plot(10*log10(alphaPower));
plot(10*log10(omegaPower));
xlabel('Channel')
ylabel('Power')
legend({'Theta','Below','Above'})
%%
figure
imagesc(10*log10(Pxx(F<40,:)'))
set(gca,'XTick',1:200:nnz(F<40),'XTickLabel',round(F(1:200:end)))
xlabel('Freq')
ylabel('channel')

%%
thetaPower_cut=thetaPower;
ff=ks.cluster_groups=='good' | ks.cluster_groups =='mua';

idx = ks.channelMap.ycoords>max(metrics.cluster_centerOfMass(ff,2));
thetaPower_cut(idx)=0;
[a,maxChan]=max(thetaPower_cut);
figure;plot(ks.channelMap.ycoords,10*log10(thetaPower(1:384)))
ylabel('LFP Power')
hold on
[a_spikes,b_spikes]=histcounts(metrics.cluster_centerOfMass(ff,2),30);

yyaxis right
plot(b_spikes(1:end-1),a_spikes)
ylabel('N Units')
%legend({'LFP Theta Power','N Units'})
xlabel('distance from tip')
%%

%%
bp = bandpass(thisDat(maxChan,:),theta_range,lfpFs);
%%
figure
plot(thisDat(maxChan,:))
hold on
plot(bp)

[pks,locs] = findpeaks(bp);
plot(locs,bp(locs),'ro')

snps=extract_snps(thisDat,locs,'win',[-400 400]);
aa_LFP=squeeze(mean(snps,3));
%%
step = 4;
mcSpike = find(ks.channelMap.ycoords>max(b_spikes),1); %only include channels with spikes
mc=385;
col = winter(ceil(mc/step));
tafig = figure('Position',[680   322   552   656]);
subplot(1,3,1)
hold on;idx =1;
tvec = [-400:400]/lfpFs;
for ii=1:step:mc
    plot(tvec,aa_LFP(ii,:)-aa_LFP(ii,1)+ii,'Color',col(idx,:));
    idx = idx+1;
end
plot(tvec,aa_LFP(maxChan,:)+maxChan,'k','LineWidth',2)
title(sprintf('max channel: %d',maxChan))
grid on
xlabel('Time')
ylabel('Distance from tip')
set(gca,'YTick',linspace(0,385,10),'YTickLabel',round(linspace(0,3840,10)))
xlim([-.2 .2]);
ylim([0 385])
%%
max_time = 400;
valid_idx = ks.spike_times>skip_seconds*30000 & ks.spike_times<(skip_seconds+max_time)*30000 & metrics.spike_is_localized; %& ismember(ks.spike_clusters,ks.clusters_good);
dd=metrics.spike_depth(valid_idx);
bins = 0:40:3840;

discrete_depth = discretize(dd,bins);
spike_times = ks.spike_times(valid_idx);
spike_times = double(spike_times);
spike_times = spike_times/30000-skip_seconds;
time_bins = 0.002; 
discrete_time = round(spike_times/time_bins)+1;
spikeMat = zeros(numel(bins),(max_time/time_bins)+1);
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
thetaPower = mean(PxxSpikes(FSpikes>theta_range(1) & FSpikes<=theta_range(2),:));


[a,maxChan_spikes]=max(thetaPower);
bp_spikes = bandpass(spikeMat(maxChan_spikes,:),theta_range,1/time_bins);
%%
figure
plot(spikeMat(maxChan_spikes,:))
hold on
plot(bp_spikes)

[pks,locs_spikes] = findpeaks(bp_spikes);
plot(locs_spikes,bp_spikes(locs_spikes),'ro')

snpsSpikes=extract_snps(spikeMat,locs_spikes,'win',[-100 100]);
aa_spikes=squeeze(mean(snpsSpikes,3));
aa_spikes=aa_spikes-mean(aa_spikes,2);
figure
imagesc(spikeMat(:,1:1000))
hold on
plot(locs_spikes(1:100),maxChan_spikes*ones(1,100),'ro')
set(gca,'YTick',linspace(0,size(spikeMat,1),10),'YTickLabel',round(linspace(min(bins),max(bins),10)))

%%
if ishandle(tafig)
figure(tafig)
else
    figure;
end
subplot(1,3,2)
hold on;idx =1;
tvec = [-100:100]/(1/time_bins);
col = winter(size(spikeMat,1));
aa_spikes = aa_spikes-aa_spikes(:,1);
for ii=1:size(spikeMat,1)
    plot(tvec,aa_spikes(ii,:)*2+ii*0.1,'Color',col(idx,:));
    idx = idx+1;
end
plot(tvec,aa_spikes(maxChan_spikes,:)*2+maxChan_spikes*0.1,'k','LineWidth',2)
title(sprintf('max depth: %d',bins(maxChan_spikes)))
grid on

set(gca,'YTick',linspace(0,size(spikeMat,1)*0.1,10),'YTickLabel',round(linspace(min(bins),max(bins),10)))
xlabel('Time')
ylabel('Distance From Tip')
ylim(size(bins)*0.1)
xlim([-.2 .2])

%% xcorr of waveforms
spike_groups = ks.spike_clusters(valid_idx);
bins = 0:80:3840;

discrete_depth = discretize(metrics.cluster_depth,bins);
avg_ACG=zeros(numel(bins),401);
n_units = zeros(1,numel(bins));
cluster_ids = ks.spike_clusters(valid_idx);

for iD=1:length(bins)
    cells = discrete_depth ==iD;
    if nnz(cells)>1
        cluids = ks.cluster_ids(cells);
        for iC=1:numel(cluids)
            this_idx = cluster_ids==cluids(iC);
            if nnz(this_idx)>100
                st_this = histcounts(spike_times(this_idx),0:time_bins:max_time);
                [xc,ll]=xcorr(st_this,200,'coeff');
                avg_ACG(iD,:)=xc+avg_ACG(iD,:);
                n_units(iD)=n_units(iD)+1;
            end
        end
    end
end
%%
avg_norm=bsxfun(@rdivide,avg_ACG,n_units');
avg_norm(:,201)=0;
avg_norm=bsxfun(@rdivide,avg_norm,sum(avg_norm,2));
avg_norm=bsxfun(@minus,avg_norm,mean(avg_norm,2));
figure
imagesc(avg_norm)
%%
if ishandle(tafig)
figure(tafig)
else
    figure;
end
subplot(1,3,3)
hold on;idx =1;
tvec = [ll]*time_bins;
col = winter(size(avg_norm,1));
for ii=1:size(avg_norm,1)
    plot(tvec,avg_norm(ii,:)*100+ii*0.1,'Color',col(idx,:));
    idx = idx+1;
end
grid on

set(gca,'YTick',linspace(0,size(avg_norm,1)*0.1,10),'YTickLabel',round(linspace(min(bins),max(bins),10)))
xlabel('Time')
ylabel('Distance From Tip')
ylim(size(bins)*0.1)
xlim([-.2 .2])
title(' Average ACG')
