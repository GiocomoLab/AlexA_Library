function [mean_waveforms,mean_template_waveforms,amplitudes,aux,clusters_good] = get_waveforms(ks_dir,imec_file,nchannels)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
chanMaps = '/oak/stanford/groups/giocomo/attialex/channelMaps';
phase1Map = fullfile(chanMaps,'phase1Map.mat');
phase0Map = fullfile(chanMaps,'phase0Map.mat');
cambridgeMap = fullfile(chanMaps,'cambridgeMap.mat');
if nchannels >300
    chanMap = phase1Map;
elseif n_channels == 277
    chanMap = phase0Map;
else
    chanMap = cambridgeMap;
end

ks = Neuropixel.KiloSortDataset(ks_dir,'channelMap',chanMap,'imecDataset',imec_file);
ks.load()
metrics = ks.computeMetrics();
nC=numel(ks.clusters_good);
clusters_good = ks.clusters_good;
mean_template_waveforms = metrics.cluster_waveform(ismember(ks.cluster_groups,'good'),:);
mean_waveforms =zeros(nC,82);
amplitudes = metrics.cluster_amplitude(ismember(ks.cluster_groups,'good'));
for iC=1:nC
    nSp = nnz(ks.spike_clusters==ks.clusters_good(iC));
    nExtract = min(200,nSp);
    snippetSet = ks.getWaveformsFromRawData('cluster_ids', ks.clusters_good(iC),'num_waveforms', nExtract, 'best_n_channels', 1, 'car', true);
    mean_waveforms(iC,:)=mean(snippetSet.data,3);
end
aux.scalefact = snippetSet.scaleToUv;
aux.timevec = snippetSet.time_ms;
end

