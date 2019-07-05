function sync_pupil_time(experiment_name,mouse_name,session_name)
if ispc()
root = 'Z:\giocomo\export\data\Projects\';
else
    root = '/oak/stanford/groups/giocomo/export/data/Projects/';
end
pupilData_name = fullfile(root,experiment_name,mouse_name,'videos',strcat(session_name,'_pupilData.mat'));
pupilTimes_name = fullfile(root,experiment_name,mouse_name,'videos',strcat(session_name,'_framedata.mat'));
vrData_name = fullfile(root,experiment_name,mouse_name,'VR',strcat(session_name,'_position.txt'));

load(pupilData_name)
load(pupilTimes_name)

%%
formatSpec = '%f%f%f%f%f%[^\n\r]';
delimiter = '\t';
fid=fopen(vrData_name);
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
vr_time = dataArray{2};
vr_time = vr_time-vr_time(1);
%%
% NIDAQ_file = 'Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I1\neuropixels_data\npI1_0415_gain_g0\npI1_0415_gain_g0_t0.nidq.bin';
% fpNIDAQ=fopen(NIDAQ_file);
% datNIDAQ=fread(fpNIDAQ,[3,Inf],'*int16');
% fclose(fpNIDAQ);
%% verify alignement, unfortnunately, the framedata.times seems to be lacking the last few frames
dpt=diff(framedata.times);
dvr=diff(vr_time(3:3:end));
figure
scatter(dpt,dvr(1:numel(dpt)),2,1:length(dpt))
title(sprintf('# diff: %i',numel(dvr)-numel(dpt)))
xlabel(session_name)
%pause
%close(gcf)
%%if good, extend pupil dat
%%
% if abs(numel(frame_times_np) - numel(frame_times_vr)) <= 1
%     idx=1:min(numel(frame_times_np),numel(frame_times_vr)); %use shorter index

%%
pupilTime = (0:0.02:max(vr_time));
vr_stamps = vr_time(3:3:end); 
%idx=1:min(numel(vr_stamps),size(pupilData,1));
idx = 1:numel(framedata.times);
%pupilData_upsampled = interp1(vr_stamps(idx),pupilData(idx,:),pupilTime);
pupilData_upsampled = interp1(framedata.times,pupilData(idx,:),pupilTime,'linear',NaN);
%figure
%plot(vr_stamps,pupilData(:,3))
%hold on
 
%plot(pupilTime,pupilData_upsampled(:,3))
%pause
%close(gcf)
save(pupilData_name,'pupilData_upsampled','pupilTime','-append');

end
