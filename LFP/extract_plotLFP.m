% myKsDir = 'F:\J1\npJ1_0525_gaincontrast_g0\npJ1_0525_gaincontrast_g0_imec0';
% channelMapFile = 'C:\code\KiloSort2\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
% im_save_dir = 'Y:\giocomo\attialex\images\LFP';
session = split(myKsDir,filesep);
session = session{end};
session = session(1:12);
im_save_dir = fullfile(im_save_dir,session);
ks = Neuropixel.KiloSortDataset(myKsDir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

nClips = 10;
clipDur = 150; % seconds

%% load data and mean subtract
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
skip_seconds = 400;
sampStarts = round(linspace(lfpFs*skip_seconds, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(lfpFs*clipDur);



mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
n=1;
thisDat = double(mmf.Data.x(:, (1:nClipSamps)+sampStarts(n)));
thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));
medianTrace = median(thisDat,1);
thisDat = thisDat-medianTrace;
%% compute psd
L = size(thisDat',1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(thisDat',[],[],[],lfpFs);
%% compute power in 3 freq bands
theta_range=[4 12];
thetaPower = mean(Pxx(F>theta_range(1) & F<=theta_range(2),:));
alphaPower = mean(Pxx(F>0.1 & F<=theta_range(1),:));
omegaPower = mean(Pxx(F>theta_range(2) & F<=20,:));
psdfig=figure();
hold on
plot(10*log10(thetaPower))
hold on
plot(10*log10(alphaPower));
plot(10*log10(omegaPower));
xlabel('Channel')
ylabel('Power')
legend({'Theta','Below','Above'})
%% plot psd across channels
powerfig=figure();
subplot(2,1,1)
imagesc(10*log10(Pxx(F<70,:)'))
set(gca,'XTick',1:200:nnz(F<70),'XTickLabel',round(F(1:200:end)))
xlabel('Freq')
ylabel('channel')
subplot(2,1,2)
mm=mean(Pxx,2);
plot(F(F<70),10*log10(mm(F<70)))

%%
thetaPower_cut=thetaPower;
ff=ks.cluster_groups=='good' | ks.cluster_groups =='mua';

%% find highest channel containing any good or mua spikes
highest_act = max(metrics.cluster_centerOfMass(ks.cluster_groups=='good',2));
idx = ks.channelMap.ycoords>highest_act;
highest_channel= find(idx,1);
thetaPower_cut(idx)=0;
[a,maxChan]=max(thetaPower_cut);

% figure;plot(ks.channelMap.ycoords,10*log10(thetaPower(1:384)))
% ylabel('LFP Power')
% hold on
% [a_spikes,b_spikes]=histcounts(metrics.cluster_centerOfMass(ff,2),30);
% 
% yyaxis right
% plot(b_spikes(1:end-1),a_spikes)
% ylabel('N Units')
% %legend({'LFP Theta Power','N Units'})
% xlabel('distance from tip')


%% bandpass channel with highest theta signal
bp = bandpass(thisDat(maxChan,:),theta_range,lfpFs);
%% find peaks in bp signal, extract raw data from all other channels (raw data)
[pks,locs] = findpeaks(bp);
tvec_lfp = [-400:1:400]/lfpFs;
snps=extract_snps(thisDat,locs,'win',[-400 400]);
%snps=extract_snps(thisDat,locs_spikesLFP,'win',[-400 400]);
aa_LFP=squeeze(mean(snps,3));
clear snps;
% figure
% plot(thisDat(maxChan,:))
% hold on
% plot(bp)
% 
% plot(locs,bp(locs),'ro')


%% plot
step = 4;
%mcSpike = find(ks.channelMap.ycoords>max(b_spikes),1); %only include channels with spikes
mc=385;
col = winter(ceil(mc/step));
tafig = figure('Position',[680   322   552   656]);

subplot(1,2,1)
hold on;idx =1;
tvec_lfp = [-400:400]/lfpFs;
chans_2_plot = 1:step:mc;
if chans_2_plot(end)<mc
    chans_2_plot = [chans_2_plot mc];
end
for ii=1:step:mc
    if ii>= highest_channel
        lc = 'k';
    else
        lc=col(idx,:);
    end
    plot(tvec_lfp,aa_LFP(ii,:)-aa_LFP(ii,1)+ii,'Color',lc);
    idx = idx+1;
end
plot(tvec_lfp,aa_LFP(maxChan,:)-aa_LFP(maxChan,1)+maxChan,'k','LineWidth',2)
title(sprintf('max channel: %d',maxChan))
grid on
xlabel('Time')
ylabel('Distance from tip')
set(gca,'YTick',linspace(0,385,10),'YTickLabel',round(linspace(0,3840,10)))
xlim([-.2 .2]);
ylim([0 385])
%%
tamapfig = figure('Position',[47         322        1185         656]);

figure(tamapfig)
subplot(1,3,1)
imagesc(flipud(aa_LFP))
set(gca,'YTick',linspace(1,385,10),'YTickLabel',round(linspace(3840,20,10)))
set(gca,'XTick',linspace(1,numel(tvec_lfp),5),'XTickLabel',linspace(min(tvec_lfp),max(tvec_lfp),5))
yline(385-highest_channel,'k');
title('LFP')
%axis image
%pcolor(aa_LFP)

%% extract spikes from raw data (for now, all detected spikes by kilosort
%% and turn it into a firing rate map split by depth
max_time = 300;
valid_idx = ks.spike_times>skip_seconds*30000 & ks.spike_times<(skip_seconds+max_time)*30000 & metrics.spike_is_localized & ismember(ks.spike_clusters,[ks.clusters_mua; ks.clusters_good]);
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
%% compute psd of this
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
locs_spikesLFP=locs_spikes*5; %turn spike time base (0.002bins) into LFP bins (*50)
locs_spikesLFP(locs_spikesLFP>(clipDur*lfpFs))=[];
% %%
% snpsLFP=extract_snps(thisDat,locs_spikesLFP,'win',[-400 400]);
% %snps=extract_snps(thisDat,locs_spikesLFP,'win',[-400 400]);
% aa_LFP=squeeze(mean(snps,3));
% figure;
% hold on
% plot(tvec_spikes,aa_spikesNorm(maxChan_spikes,:))
% 
% [ff,minloc]=min(abs(bins(maxChan_spikes)-ks.channelMap.ycoords));
% plot(tvec_lfp,aa_LFP(minloc,:))

%% plot a small excerpt

% figure
% plot(spikeMat(maxChan_spikes,:))
% hold on
% plot(bp_spikes)
% 
% plot(locs_spikes,bp_spikes(locs_spikes),'ro')

snpsSpikes=extract_snps(spikeMat,locs_spikes,'win',[-100 100]);
tvec_spikes = [-100:100]*time_bins;
aa_spikes=squeeze(mean(snpsSpikes,3));
clear snpsSpikes;
aa_spikesNorm=bsxfun(@rdivide,aa_spikes,mean(aa_spikes,2));
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


hold on
[~,max_loc]=max(aa_spikes,[],2);
hold on
plot((max_loc(end:-1:1)),1:numel(max_loc),'ro')
%%
snpsLFP=extract_snps(thisDat,locs_spikesLFP(1:100),'win',[-400 400]);
%snps=extract_snps(thisDat,locs_spikesLFP,'win',[-400 400]);
aa_LFP=squeeze(mean(snpsLFP,3));
clear('snpsLFP')
%%
spikelfpfig=figure('Position',[1436         161         392         817]);
chan_idx=[-20:5:20]+maxChan_spikes;
chan_idx(chan_idx<1)=[];
cntr=0;
for ii=numel(chan_idx):-1:1
    cntr = cntr+1;
    try
subplot(numel(chan_idx),1,cntr)
hold on
plot(tvec_spikes,aa_spikes(chan_idx(ii),:))
title(bins(chan_idx(ii))-bins(maxChan_spikes))

[ff,minloc]=min(abs(bins(chan_idx(ii))-ks.channelMap.ycoords));
yyaxis right
plot(tvec_lfp,aa_LFP(minloc,:))
grid on
axis tight
    catch
    end
end
legend({'spikes','LFP'})


%%
figure(tamapfig)
subplot(1,3,2)
%aa_spikes_norm = 
imagesc(flipud(aa_spikesNorm),[.05 8])
set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
%yline(385-highest_channel,'k');
xline(numel(tvec_spikes)/2+.5)
title('triggered spikes')

%%
if ishandle(tafig)
figure(tafig)
else
    figure;
end
subplot(1,2,2)
hold on;idx =1;
tvec = [-100:100]/(1/time_bins);
col = winter(size(spikeMat,1));
aa_spikes_plot = aa_spikes-aa_spikes(:,1);
for ii=1:size(spikeMat,1)
    plot(tvec,aa_spikes_plot(ii,:)*2+ii*0.1,'Color',col(idx,:));
    idx = idx+1;
end
plot(tvec,aa_spikes_plot(maxChan_spikes,:)*2+maxChan_spikes*0.1,'k','LineWidth',2)
title(sprintf('max depth: %d',bins(maxChan_spikes)))
grid on

set(gca,'YTick',linspace(0,size(spikeMat,1)*0.1,10),'YTickLabel',round(linspace(min(bins),max(bins),10)))
xlabel('Time')
ylabel('Distance From Tip')
ylim(size(bins)*0.1)
xlim([-.2 .2])
%%

%% xcorr of spike times
max_time = 1800;
valid_idx = ks.spike_times>skip_seconds*30000 & ks.spike_times<(skip_seconds+max_time)*30000 & metrics.spike_is_localized & ismember(ks.spike_clusters,[ks.clusters_good]);
spike_times = ks.spike_times(valid_idx);
spike_times = double(spike_times);
spike_times = spike_times/30000-skip_seconds;
time_bins = 0.002; 
spike_groups = ks.spike_clusters(valid_idx);
bins_xcorr = 0:100:3840;

discrete_depth = discretize(metrics.cluster_depth,bins_xcorr);
avg_ACG=zeros(numel(bins_xcorr),401);
n_units = zeros(1,numel(bins_xcorr));
cluster_ids = ks.spike_clusters(valid_idx);

for iDepth=1:length(bins_xcorr)
    cells = discrete_depth ==iDepth;
    if nnz(cells)>1
        cluids = ks.cluster_ids(cells);
        for iC=1:numel(cluids)
            this_idx = cluster_ids==cluids(iC);
            if nnz(this_idx)>100
                st_this = histcounts(spike_times(this_idx),0:time_bins:max_time);
                [xc,ll]=xcorr(st_this,200,'coeff');
                avg_ACG(iDepth,:)=xc+avg_ACG(iDepth,:);
                n_units(iDepth)=n_units(iDepth)+1;
            else
                %sprintf('few spikes: %d %d',nnz(this_idx),cluids(iC))
            end
        end
    end
end
%%
tvec = ll*time_bins;

avg_norm=bsxfun(@rdivide,avg_ACG,n_units');
avg_norm(:,201)=0;
%avg_norm=bsxfun(@rdivide,avg_norm,sum(avg_norm,2));
%avg_norm=bsxfun(@minus,avg_norm,mean(avg_norm,2));
figure(tamapfig)
subplot(1,3,3)
%aa_spikes_norm = 
imagesc(flipud(avg_norm))
set(gca,'YTick',linspace(1,size(avg_ACG,1),10),'YTickLabel',round(linspace(max(bins_xcorr),min(bins_xcorr),10)))
set(gca,'XTick',linspace(1,numel(ll),5),'XTickLabel',linspace(min(tvec),max(tvec),5))
yline(385-highest_channel,'k');
title('triggered spikes')
%% save figures
suffix='_CAR';
saveas(spikefig,strcat(im_save_dir,'spikes',suffix,'.png'));
close(spikefig)
saveas(tamapfig,strcat(im_save_dir,'heatmaps',suffix,'.png'));
close(tamapfig)
saveas(tafig,strcat(im_save_dir,'LFP',suffix,'.png'));
close(tafig)
saveas(powerfig,strcat(im_save_dir,'Power',suffix,'.png'));
close(powerfig)
saveas(spikelfpfig,strcat(im_save_dir,'spikeLFP',suffix,'.png'));
close(spikelfpfig);
saveas(psdfig,strcat(im_save_dir,'LFPPower',suffix,'.png'));
close(psdfig)
%%
% if ishandle(tafig)
% figure(tafig)
% else
%     figure;
% end
% subplot(1,3,3)
% hold on;idx =1;
% tvec = [ll]*time_bins;
% col = winter(size(avg_norm,1));
% for ii=1:size(avg_norm,1)
%     plot(tvec,avg_norm(ii,:)*100+ii*0.1,'Color',col(idx,:));
%     idx = idx+1;
% end
% grid on
% 
% set(gca,'YTick',linspace(0,size(avg_norm,1)*0.1,10),'YTickLabel',round(linspace(min(bins),max(bins),10)))
% xlabel('Time')
% ylabel('Distance From Tip')
% ylim(size(bins)*0.1)
% xlim([-.2 .2])
% title(' Average ACG')
