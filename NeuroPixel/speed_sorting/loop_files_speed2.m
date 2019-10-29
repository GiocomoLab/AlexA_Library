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
    if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
        session_name{end+1}=sn(iS).name(1:end-4);
    end
end
% all the values we ever use
gains_all = [0.8 0.7 0.6 0.5 0.2];
contrasts_all = [100 50 20 10 5 2 0];
aggregateData = struct();

%% iterate over sessions
parpool(8)
savedir = fullfile(OAK,'attialex','xcorrFlat');
if ~isfolder(savedir)
    mkdir(savedir)
end
parfor session_num = 1:numel(session_name)  
     
    % load data
    fprintf('session %d/%d: %s\n',session_num,numel(session_name),session_name{session_num});

try
    fit_speed2(fullfile(data_dir,session_name{session_num}),savedir,params)
catch
    sprintf('failed for %s',session_name{session_num})
end
end

