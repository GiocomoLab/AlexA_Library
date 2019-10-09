% paths
%addpath(genpath('/home/users/attialex/neuropixel-utils'))
%channelMapFile = fullfile('/home/users/attialex/neuropixel-utils/map_files','neuropixPhase3B1_kilosortChanMap.mat');
addpath(genpath('F:\code\neuropixel-utils'))
channelMapFile = fullfile('F:\code','neuropixel-utils','map_files','neuropixPhase3B1_kilosortChanMap.mat');

%dataset to load
% ksdirs=dir('Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\**\*imec0.ap.bin');
% session = split(myKsDir,filesep);
% session = session{end};
% session = session(1:12);
myKsDir='Z:\giocomo\export\data\Projects\AlexA_NP\AA_190906_7\neuropixels_data\AA_190906_7_1004_gaincontrast10_1_g0\AA_190906_7_1004_gaincontrast10_1_g0_imec0';
%%
ks = Neuropixel.KiloSortDataset(myKsDir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();
%% all
metrics.plotClusterWaveformAtCenterOfMass()
metrics.plotDriftmap('driftThreshold', 4, 'spikeAmpQuantile', 0.8);


cluster_ids = metrics.cluster_ids(1:5:metrics.nClusters);
metrics.plotDriftmap('tsi', tsi, 'exciseRegionsOutsideTrials', true, 'cluster_ids', cluster_ids);

%% cluster specific stuff
for iC=1:clusters
     [CGR,b]=CCG(sp.st,cluID,'binSize',[0.001],'duration',[0.2]);

    metrics.plotClusterImage(cluster_id, 'best_n_channels', 20);
snippetSet = ks.getWaveformsFromRawData('cluster_ids', ks.clusters_good(1), ...
    'num_waveforms', 500, 'best_n_channels', 20, 'car', true, ...
    'subtractOtherClusters', true);
end