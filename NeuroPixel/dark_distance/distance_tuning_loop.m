root = 'F:\attialex\';
root = '/oak/stanford/groups/giocomo/attialex/'
files = dir(fullfile(root,'NP_DATA_corrected','np*dark*'));
savepath = fullfile(root,'distance_tuning');

if ~isfolder(savepath)
    mkdir(savepath)
end

ops = load_default_opt;

ops.num_shuf = 400;
ops.dark = true;
ops.SpatialBin = 5;
ops.max_lag = 800;

%%

for iF=1:numel(files)
    save_name = fullfile(savepath,files(iF).name);
    if isfile(save_name)
        disp('exists')
        continue
    end
    data = load(fullfile(files(iF).folder,files(iF).name));
    
    data_out = calc_distance_tuning(data,data.sp.cids(data.sp.cgs==2),ops);
    
    save_name = fullfile(savepath,files(iF).name);
    save(save_name,'data_out')
end