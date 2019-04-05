addpath(genpath('C:\code\AlexA_Library'))
root=fullfile('Y:','giocomo','attialex','NP_DATA');
Files= dir(fullfile(root,'*.mat'));

for iF=1:length(Files)
    try
        clearvars -except iF root Files
    load(fullfile(root,Files(iF).name));
    pairwise_correlation
    catch ME
        display(['did not work for ' Files(iF).name])
        display(ME.Identifier)
    end
    drawnow
end