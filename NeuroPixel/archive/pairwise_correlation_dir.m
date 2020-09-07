addpath(genpath('C:\code\AlexA_Library'))
root=fullfile('Z:','giocomo','attialex','NP_DATA');
Files= dir(fullfile(root,'npI5_0417_baseline_1.mat'));

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