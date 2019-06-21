myKsDirExternal ='Y:\giocomo\export\data\Projects\RandomForage_NPandH3\Rome\Rome_0423_2_g0\0423_2_g0_imec0';
myKsDirInternal = 'Y:\giocomo\export\data\Projects\RandomForage_NPandH3\Rome\Rome_0422_g0\0422_g0_imec0'; %internal
myKsDir = myKsDirExternal;
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

%chanMap = readNPY(fullfile(myKsDir, 'channel_map.npy'));
chanMap = 1:383;
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

% plot LFP power
dispRange = [0 100]; % Hz
marginalChans = [10:50:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};
figure;
plotLFPpower(F, allPowerEst, dispRange, 100:105, freqBands,true(1,nChansInFile));
%%

nClips = 10;
clipDur = 200; % seconds

d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
skip_seconds = 100;
sampStarts = round(linspace(lfpFs*skip_seconds, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(lfpFs*clipDur);



mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
n=1;
thisDat = double(mmf.Data.x(:, (1:nClipSamps)+sampStarts(n)));
thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));

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

figure
img = 10*log10(Pxx(F<40,:)');
imagesc(flipud(img))
set(gca,'XTick',1:200:nnz(F<40),'XTickLabel',round(F(1:200:end)))
set(gca,'YTick',0:50:350,'YTickLabel',383:-50:0)
xlabel('Freq')
ylabel('channel')

thetaPower_cut=thetaPower;

[a,maxChan]=max(thetaPower_cut);



%% bandpass channel with highest theta signal
bp = bandpass(thisDat(maxChan,:),theta_range,lfpFs);
% find peaks in bp signal, extract raw data from all other channels (raw data)
[pks,locs] = findpeaks(bp);
snps=extract_snps(thisDat,locs,'win',[-400 400]);
aa_LFP=squeeze(mean(snps,3));



%% plot
step = 4;
mc=385;
col = winter(ceil(mc/step));
tafig = figure('Position',[680   322   552   656]);


hold on;idx =1;
tvec = [-400:400]/lfpFs;
chans_2_plot = 1:step:mc;
if chans_2_plot(end)<mc
    chans_2_plot = [chans_2_plot mc];
end
for ii=1:step:mc
    
        lc=col(idx,:);
    
    plot(tvec,aa_LFP(ii,:)-aa_LFP(ii,1)+ii,'Color',lc);
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