% paths
%addpath(genpath('/home/users/attialex/neuropixel-utils'))
%channelMapFile = fullfile('/home/users/attialex/neuropixel-utils/map_files','neuropixPhase3B1_kilosortChanMap.mat');
%addpath(genpath('C:\code\neuropixel-utils'))
%channelMapFile = fullfile('C:\code','neuropixel-utils','map_files','neuropixPhase3B1_kilosortChanMap.mat');
addpath(genpath('/home/users/attialex/neuropixel-utils'))
channelMapFile = fullfile('/home/users/attialex/neuropixel-utils/map_files','neuropixPhase3B1_kilosortChanMap.mat');
%ksdirs = dir(fullfile('/oak/stanford/groups/giocomo','export','data','Projects','RandomForage_NPandH3','ProcessedData','/Hanover*/*ap.bin'));
%ksdirs = dir(fullfile('/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/**','*ap.bin'))
%ksdirs=dir('Y:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\**\*imec0.ap.bin');
ksdirs = dir(fullfile('/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/AA_*/neuropixels_data/AA*/AA*/*.ap.bin'));
%dataset to load
% ksdirs=dir('Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\**\*imec0.ap.bin');
% session = split(myKsDir,filesep);
% session = session{end};
% session = session(1:12);
%ksDirs = dir('F:\Alex\**\*imec0.ap.bin');
im_root = fullfile('/oak/stanford/groups/giocomo','images');
%im_root = 'F:\images';

%%
for iDir = 1:numel(ksDirs)
%myKsDir='F:\Alex\AA_190906_7\neuropixels_data\AA_190906_7_1004_gaincontrast10_1_g0\AA_190906_7_1004_gaincontrast10_1_g0_imec0';
myKsDir = ksDirs(iDir).folder;
%%
session = split(myKsDir,filesep);
session = session{end};
%session = session(1:12);
im_save_dir = fullfile(im_root,session);
if ~isfolder(im_save_dir)
    mkdir(im_save_dir)
end

ks = Neuropixel.KiloSortDataset(myKsDir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();

centermass = figure('Position',[ 680   112   381   866]);
metrics.plotClusterWaveformAtCenterOfMass('cluster_IDs',ks.clusters_good)
driftmap_1 = figure('Position',[680   558   925   420]);
metrics.plotDriftmap('driftThreshold', 4, 'spikeAmpQuantile', 0.8);


cluster_ids = ks.clusters_good(1:3:end);
driftmap_2 = figure('Position',[680   558   925   420]);

metrics.plotClusterDriftmap('cluster_IDs', cluster_ids);
driftmap_3 = figure('Position',[680   558   925   420]);
metrics.plotClusterDriftmap('cluster_IDs', cluster_ids,'showSmooth', true, 'showIndividual', false, 'smoothWidthSeconds', 50)
saveas(centermass,fullfile(im_save_dir,'centermass.png'))
close(centermass)
for ii=1:3
    eval(sprintf('cf=driftmap_%d;',ii));
    saveas(cf,fullfile(im_save_dir,sprintf('driftmap_%d.png',ii)))
    close(cf)
end
% 
% %% cluster specific stuff
n_good = numel(ks.clusters_good);
figure('Position',[378   278   862   700])
gc=ismember(ks.cluster_groups,'good');
amp = metrics.cluster_amplitude(gc);
[~,sid]=sort(amp,'descend');
for iC=1:10
    cluid = ks.clusters_good(sid(iC));
    idx=ks.spike_clusters==cluid;
    bb=double(ks.spike_times(idx))/ks.sample_rate;

    subplot(3,2,1)
    plot(bb,metrics.spike_amplitude(idx),'.')
    
%     [CGR,b]=CCG(bb,ones(nnz(idx),1));
%     subplot(3,1,2)
%     bar(b,CGR,'EdgeColor','none')
    subplot(3,2,2)
    [CGR,b]=CCG(bb,ones(nnz(idx),1),'binSize',0.001,'duration',.100);

    bar(b,CGR,'EdgeColor','none')
    
    subplot(3,2,[3 5])
    metrics.plotClusterImage(cluid, 'best_n_channels', 20);
    pbaspect([1.5,4,1])

    snippetSet = ks.getWaveformsFromRawData('cluster_ids', cluid, ...
    'num_waveforms', 100, 'best_n_channels', 20, 'car', true, ...
    'subtractOtherClusters', false);
    subplot(3,2,[4 6])
    snippetSet.plotAtProbeLocations('alpha',0.8,'showChannelLabels',false)
    pbaspect([1.5,4,1])
    %saveas(gcf,fullfile(im_save_dir,sprintf('%d.pdf',cluid)))
    gcf
    print(fullfile(im_save_dir,sprintf('%d.png',cluid)),'-dpng','-r300')
end
close(gcf)
end