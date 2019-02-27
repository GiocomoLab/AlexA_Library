%function create_h5_files(DATADIR)
addpath('~/AlexA_Library/suite2p_helpers')
addpath('~/AlexA_Library/suite2p_helpers/rigid_registration')


DATADIR = '/oak/stanford/groups/giocomo/attialex/TEST';


content = dir(DATADIR);

for iF = 1:length(content)
    if content(iF).isdir
        rigid_reg_for_dir(fullfile(DATADIR,content(iF).name))
    end
end
