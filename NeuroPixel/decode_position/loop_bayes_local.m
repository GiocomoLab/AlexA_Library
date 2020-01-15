addpath(genpath('C:\code\AlexA_Library'));
params = readtable('UniversalParams.xlsx');

% where to find data and save images
%data_dir = 'Z:\giocomo\attialex\NP_DATA\';\

data_dir='Z:\giocomo\attialex\NP_DATA';
%session_name = {'AA5_190809_gain_1'};
session_name = {};
sn = dir(fullfile(data_dir,'np*_gain*.mat'));
for iS = 1:numel(sn)
    if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
        session_name{end+1}=sn(iS).name(1:end-4);
    end
end
% all the values we ever use
load('Z:\giocomo\attialex\cellInfoGain.mat')
%% iterate over sessions
parpool(4)
savedir = fullfile('Z:\giocomo\attialex\','decodeBayesv3');

if ~isfolder(savedir)
    mkdir(savedir)
end

parfor session_num = 1:numel(session_name)  
     
    % load data
    fprintf('session %d/%d: %s\n',session_num,numel(session_name),session_name{session_num});

try
    decode_bayes(fullfile(data_dir,session_name{session_num}),savedir,cell_info);
catch
    sprintf('failed for %s',session_name{session_num})
end
end

