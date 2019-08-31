root = 'Z:\giocomo\attialex\NP_DATA';
Files = dir(fullfile(root,'*_gain_*.mat'));
MERGED=struct;
for iF=1:numel(Files)
    clearvars -except root Files iF MERGED
    fn = fullfile(root,Files(iF).name);
    dataset = load(fn);
    similarity_gain
    MERGED(iF).stability = stability;
    MERGED(iF).similarity = similarity;
    MERGED(iF).contrast = dataset.trial_contrast;
    MERGED(iF).gain = dataset.trial_gain;
    MERGED(iF).avgMaps = avgMaps;
    MERGED(iF).avgMapsSimilarity = avgMapSimilarity;
    MERGED(iF).maps_smoothed = maps_smoothed;
end

%%
