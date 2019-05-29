root = 'Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\J1\videos';
parts = split(root,'\');
exp=parts{end-2};
animal=parts{end-1};
fn=dir([root filesep '*.avi']);

for iF =1%:length(fn)
    if ~contains(fn(iF).name,'dark')
    analyzeFaceVideo(fullfile(root,fn(iF).name));
    sync_pupil_time(exp,animal,fn(iF).name(1:end-4));
    end
end