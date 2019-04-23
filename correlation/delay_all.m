root = 'Z:\giocomo\attialex\NP_DATA\';
addpath(genpath('C:\code\spikes\analysis'),'-end');
matfiles = dir([root '*.mat']);

all_depths={};
all_delays={};
all_CGR={};
for iF=1:length(matfiles)
    clear sp
    load(fullfile(root,matfiles(iF).name))
    try
    delay_along_probe
    set(gcf,'Name',matfiles(iF).name)
    saveas(gcf,fullfile('C:','temp',[matfiles(iF).name '_delay.png']))
    close(gcf)
    all_depths{end+1}=mean_depth;
    all_delays{end+1}=delays;
    all_CGR{end+1}=CGR;
    catch ME
        sprintf('error with file %s \n',matfiles(iF).name);
        disp(ME.identifier)
    end
end
%%