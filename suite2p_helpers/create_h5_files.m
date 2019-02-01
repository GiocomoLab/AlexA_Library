%function create_h5_files(DATADIR)
addpath('~/AlexA_Library/suite2p_helpers')

DATADIR = '/oak/groups/giocomo/attialex/DATA/TEST';


content = dir(DATADIR);

for iF = 1:length(content)
    if content(iF).isDir
        h5_for_dir(fullfile(DATADIR,content(iF).name))
    end
end
