pause(1800);

addpath(genpath('C:\code\KiloSort2')) % path to kilosort folder
addpath('C:\code\npy-matlab')
addpath(genpath('C:\code\AlexA_Library'))
addpath('C:\code\SiliconProbeCode\')

directories={'F:\J1\npJ1_0525_gaincontrast_g0\npJ1_0525_gaincontrast_g0_imec0'};
for dd = 1:length(directories)
    tic
    try
    run_ks2(directories{dd});
    catch ME
        fprintf('failed for %s \n',directories{dd})
        fprintf('Error: %s \n',ME.identifier)
    end
    toc
end