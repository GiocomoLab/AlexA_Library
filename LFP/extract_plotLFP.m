myKsDir = 'F:\J1\npJ1_0520_baseline_g0\npJ1_0520_baseline_g0_imec0';
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

nClips = 10;
clipDur = 50; % seconds

% load nClips one-sec samples
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
sampStarts = round(linspace(lfpFs*10, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(lfpFs*clipDur);



mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
n=1;
thisDat = double(mmf.Data.x(:, (1:nClipSamps)+sampStarts(n)));
thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));
%%
L = size(thisDat',1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(thisDat',[],[],NFFT,lfpFs);
%%

thetaPower = mean(Pxx(F>3 & F<13,:));
%%
[a,maxChan]=max(thetaPower);
figure;plot(ks.channelMap.ycoords,thetaPower(1:384))
hold on
ff=ks.cluster_groups=='good' | ks.cluster_groups =='mua';
[a,b]=histcounts(metrics.cluster_centerOfMass(ff,2),30);
yyaxis right
plot(b(1:end-1),a)
legend({'LFP Theta Power','N Units'})
xlabel('distance from tip')
%%

%%
bp = bandpass(thisDat(maxChan,:),[8 13],2500);

figure
plot(thisDat(maxChan,:))
hold on
plot(bp)

[pks,locs] = findpeaks(bp);
plot(locs,bp(locs),'ro')

snps=extract_snps(thisDat,locs,'win',[-200 200]);
aa=squeeze(mean(snps(:,:,1:10),3));

step = 4;
col = parula(ceil(nChansInFile/step));
figure;hold on;idx =1;for ii=1:step:nChansInFile;plot(aa(ii,:),'Color',col(idx,:));idx = idx+1;end;