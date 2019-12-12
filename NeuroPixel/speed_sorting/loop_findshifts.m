

ops.factors = -.25:0.01:.25;
ops.edges = 0:2:400;
%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = 3:20;
ops.TimeBin = 0.02;
fi = gausswin(5);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = false;
OAK='/oak/stanford/groups/giocomo/';
%% savedir = 
savedir = fullfile(OAK,'attialex','speed_shift');
if ~isfolder(savedir)
    mkdir(savedir);
end
mf = matfile(fullfile(savedir,'parameters'));
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


for iF=1:numel(filenames)
    data = load(fullfile(data_dir,filenames{iF}));
    data_out = findBestShifts(data,ops);
    [~,filenames,~]=fileparts(filenames{iF});
    mf =matfile(fullfile(savepath,filenames));
    mf.stability = data_out.all_stability;
    mf.region = data_out.region;
end
