

ops.factors = -.25:0.01:.25;
ops.edges = 0:2:400;
%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = 3:20;
ops.TimeBin = 0.02;
ops.idx = [10:390];
fi = gausswin(5);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = false;
OAK='/oak/stanford/groups/giocomo/';
%% savedir = 
savedir = fullfile(OAK,'attialex','speed_shiftFilter');
if ~isfolder(savedir)
    mkdir(savedir);
end
mf = matfile(fullfile(savedir,'parameters'),'Writable',true);
mf.ops = ops;


%% find files

data_dir=fullfile(OAK,'attialex','NP_DATA');
%session_name = {'AA5_190809_gain_1'};
filenames = {};
sn = dir(fullfile(data_dir,'*.mat'));
for iS = 1:numel(sn)
    if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
        filenames{end+1}=sn(iS).name(1:end-4);
    end
end

%%
p=gcp('nocreate');
if isempty(p)
    parpool(12);
end

%%


parfor iF=1:numel(filenames)
    try
    data = load(fullfile(data_dir,filenames{iF}));
    ops_here = ops;
    %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
    data_out = findBestShifts(data,ops);
    [~,session_name,~]=fileparts(filenames{iF});
    mf =matfile(fullfile(savedir,session_name));
    mf.stability = data_out.all_stability;
    mf.region = data_out.region;
    catch ME
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
    end
end
