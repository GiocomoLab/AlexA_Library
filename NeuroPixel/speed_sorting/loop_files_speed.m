addpath(genpath('/home/users/attialex/AlexA_Library'));
params = readtable('UniversalParams.xlsx');
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
%data_dir = 'Z:\giocomo\attialex\NP_DATA\';\
run('/home/users/attialex/AlexA_Library/default_paths.m')
data_dir=fullfile(OAK,'attialex','NP_DATA');
%session_name = {'AA5_190809_gain_1'};
session_name = {};
sn = dir(fullfile(data_dir,'*_gain_*.mat'));
for iS = 1:numel(sn)
    if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback'))
        session_name{end+1}=sn(iS).name(1:end-4);
    end
end
% all the values we ever use
gains_all = [0.8 0.7 0.6 0.5 0.2];
contrasts_all = [100 50 20 10 5 2 0];
aggregateData = struct();

%% iterate over sessions
for session_num = 1:numel(session_name)  
        save_figs=true;
        image_save_dir = strcat('/oak/stanford/groups/giocomo/attialex/Images/',...
    session_name{session_num},'speed_sort');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end

    % load data
    fprintf('session %d/%d: %s\n',session_num,numel(session_name),session_name{session_num});
    clearvars -except sess* aggre* data* save_figs image_save_dir
    load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));
    shift_vs_speed
    aggregateData(session_num).models = models;
    aggregateData(session_num).session = session_name{session_num};
    aggregateData(session_num).region = region(sp.cgs==2);
    aggregateData(session_num).c_coeff = c_coeff;
    aggregateData(session_num).p_vals = p_vals;
    aggregateData(session_num).clu_reg = clu_reg;
end

save('/oak/stanford/groups/giocomo/attialex/shift_vs_speed.mat',aggregateData);