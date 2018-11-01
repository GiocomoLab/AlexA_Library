%myKsDir='F:\F3\1018_contrasttrack_gainchanges_1';
%sp = loadKSdir(myKsDir);
%%
gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(2).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 277;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st(sp.clu==2)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==2);

wf = getWaveForms(gwfparams);

figure; 
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;

[a,b]=max(wf.waveFormsMean,[],2);
[ii,jj]=max(a);
figure
plot(squeeze(wf.waveFormsMean(1,b(jj),:)))