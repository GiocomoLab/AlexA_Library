root = 'Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I1\videos';

fn=dir([root filesep '*.avi']);

for iF =25:length(fn)
    if ~contains(fn(iF).name,'dark')
    analyzeFaceVideo(fullfile(root,fn(iF).name));
    end
end