%check if this experiment already in DB

%if not, add experiment to db

%get prefix for cell #
ks = Neuropixel.KiloSortDataset(myKsDir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();

for iC=1:clusters

snippetSet = ks.getWaveformsFromRawData('cluster_ids', ks.clusters_good(1), ...
    'num_waveforms', 500, 'best_n_channels', 20, 'car', true, ...
    'subtractOtherClusters', true);

ks = Neuropixel.KiloSortDataset(myKsDir,'channelMap',channelMapFile);
ks.load()
metrics = ks.computeMetrics();

%add cluster info to db

% save mean waveform data cell_prefix cluID (4 digits, 4 digits) 10000 0000
% and channel coordinate data


%%