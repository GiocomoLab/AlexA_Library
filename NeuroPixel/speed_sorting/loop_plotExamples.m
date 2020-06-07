shift_path = '/Volumes/Samsung_T5/attialex/speed_filtered_correctedData';
data_path = '/Volumes/Samsung_T5/attialex/NP_Data_corrected';
ops = load(fullfile(shift_path,'parameters.mat'));
ops = ops.ops;
ops.speedWindow = [-10 -1]; % in cm
ops.midpoints = ops.edges(1:end-1)*.5+ops.edges(2:end)*.5;






%%
gain = 0.8;
contrast = 100;
regions = {'VISp'};
filenames = {};
triggers = {};
ops.savepath =fullfile('/Volumes/Samsung_T5/attialex','shift_images_corrected',regions{1});
if ~isfolder(ops.savepath)
    mkdir(ops.savepath)
end

for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,data_path);
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
%%

for iF=1:numel(filenames)
    [~,sn]=fileparts(filenames{iF});
    if ~isfile(fullfile(shift_path,[sn,'.mat']))
        disp('no shift')
        continue
    end
    data = load(filenames{iF});
    shift_data = load(fullfile(shift_path,sn));
    trials = -15:13;
    for iT=1%:numel(triggers{iF})
        ops.trials = triggers{iF}(iT)+trials;
        
        
        
        ops.bl_pre = 1:15;
        ops.gain_trials = 16:19;
        ops.bl_post = 20:29;
        plotExampleCells(data,shift_data,regions{1},sn,ops);
    end
end
