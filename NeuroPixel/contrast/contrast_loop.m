root = 'Z:\giocomo\attialex\NP_DATA';
Files = dir(fullfile(root,'*_contrast_*.mat'));
MERGED=struct;
for iF=1:15
    clearvars -except root Files iF MERGED
    fn = fullfile(root,Files(iF).name);
    contrast = load(fn);
    similarity_contrast
    MERGED(iF).stability = stability;
    MERGED(iF).similarity = similarity;
    MERGED(iF).contrast = contrast.trial_contrast;
    MERGED(iF).gain = contrast.trial_gain;
    MERGED(iF).avgMaps = avgMaps;
    MERGED(iF).avgMapsSimilarity = avgMapSimilarity;
    MERGED(iF).maps_smoothed = maps_smoothed;
    MERGED(iF).gain_response = gain_response;
    MERGED(iF).gain_stability = gain_stability;
    MERGED(iF).xcorr_pre_post = xcorr_pre_post;
    MERGED(iF).xcorr_pre_gain = xcorr_pre_gain;
end

%%
