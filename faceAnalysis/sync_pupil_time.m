function sync_pupil_time(experiment_name,mouse_name,session_name)
root = 'Z:\giocomo\export\data\Projects\';
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

%%
% NIDAQ_file = 'Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I1\neuropixels_data\npI1_0415_gain_g0\npI1_0415_gain_g0_t0.nidq.bin';
% fpNIDAQ=fopen(NIDAQ_file);
% datNIDAQ=fread(fpNIDAQ,[3,Inf],'*int16');
% fclose(fpNIDAQ);
%% verify alignement, unfortnunately, the framedata.times seems to be lacking the last few frames
dpt=diff(framedata.times);
dvr=diff(dataArray{2}(3:3:end));
figure
plot(dpt,dvr(1:numel(dpt)),'.')
pause
close(gcf)
%%if good, extend pupil dat
%%
% if abs(numel(frame_times_np) - numel(frame_times_vr)) <= 1
%     idx=1:min(numel(frame_times_np),numel(frame_times_vr)); %use shorter index

%%
pupilTime = (0:0.02:max(dataArray{2}));
vr_stamps = dataArray{2}(3:3:end);
pupilData_upsampled = interp1(vr_stamps,pupilData,pupilTime);

figure
plot(dataArray{2}(3:3:end),pupilData(:,3))
hold on
 
plot((0:0.02:max(dataArray{2})),pupilData_upsampled(:,3))
pause
close(gcf)
save(pupilData_name,'pupilData_upsampled','pupilTime','-append');

end
