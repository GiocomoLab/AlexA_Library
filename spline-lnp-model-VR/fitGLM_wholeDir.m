%addpath(genpath('F:\code\cortexlab_spikes\analysis\'))

%matfiles = dir('/Users/attialex/NP_DATA_2/*mismatch*.mat');
OAK = '/Volumes/T7/attialex';
OAK='/oak/stanford/groups/giocomo/attialex';

matfiles = dir(fullfile(OAK,'NP_DATA_corrected/*.mat'));
savedir = fullfile(OAK,'glmFits');
if ~isfolder(savedir)
    mkdir(savedir)
end

%%

parfor iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    if isfile(fullfile(savedir,matfiles(iF).name))
        continue
    end
    try
    good_cells=data_out.sp.cids(data_out.sp.cgs==2);
    trials = 1:(max(data_out.trial));
    trials = trials(data_out.trial_contrast == 100 & data_out.trial_gain==1);
    glmData = fitGLM(data_out,trials,good_cells);
    mf = matfile(fullfile(savedir,matfiles(iF).name))
    mf.glmData = glmData;
    %save(fullfile(savedir,matfiles(iF).name),'glmData');
    catch ME
        disp(ME.message)
    end
end
    