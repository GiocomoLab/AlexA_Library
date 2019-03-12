root=fullfile('/oak','stanford','groups','giocomo','attialex','NP_DATA');
Files= dir(fullfile(root,'*.mat'));
addpath('..')
addpath('../NeuroPixel')
params = readtable('UniversalParams.xlsx');
close all;

for iF=1:length(Files)
    mkdir(root,'Plots');
    [~,session,~]=fileparts(Files(iF).name);
    mkdir(fullfile(root,'Plots'),session);
    plot_path=fullfile(root,'Plots',session);
    try
        fprintf('Now working on %s \n',Files(iF).name )
        load(fullfile(root,Files(iF).name));
        glm_playground
        fprintf('Done for %s, now saving. \n',Files(iF).name)
        save(fullfile(root,Files(iF).name),'glmData','-append')
    catch ME
        fprintf('Failed for %s \n',Files(iF).name)
        warning(ME.message)
    end
end
