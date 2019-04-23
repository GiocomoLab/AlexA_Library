function [vr_data2p] = sync_2p_vr(logfile, infostruct)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
first_frame_2p = min(infostruct.frame);
last_frame_2p = max(infostruct.frame);
fps_2p=infostruct.resfreq*2/infostruct.sz(1);

n_frames_with_stim = last_frame_2p-first_frame_2p;
valid_frames_2p=false(1,last_frame_2p);
valid_frames_2p(first_frame_2p:last_frame_2p)=true;



delimiter = ',';
fid = fopen(logfile,'r');
header = fgetl(fid);
headers = strtrim(split(header,','));
formatSpec='%s%s';
for ii = 1:length(headers)-2
    formatSpec=[formatSpec '%f'];
end
        
%formatSpec = '%s%s%f%f%f%[^\n\r]';
formatSpec =[formatSpec '%[^\n\r]'];
%%
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fid);
%%
vr_time=zeros(1,size(dataArray{2},1));
nvr_frames=size(dataArray{2},1);
for ii=1:nvr_frames
    vr_time(ii)=posixtime(datetime(dataArray{2}{ii},'InputFormat','HH:mm:ss.S'));
end
%%
vr_time=vr_time-vr_time(1);
vr_data = struct();
for ii=3:length(headers)
vr_data.(headers{ii}) = cat(1,dataArray{ii});
end

%%
time_vec_2p=zeros(1,n_frames_with_stim);
pulse=1;
for ii=1:length(time_vec_2p)
    current_frame = first_frame_2p+(ii-1);
    n_pulse_frame = nnz(infostruct.frame==current_frame);
    start= pulse;
    stop = pulse+(n_pulse_frame-1);
    pulse = pulse+n_pulse_frame;
    time_vec_2p(ii)=mean(vr_time(start:stop));
end
%%
vr_data2p=struct();
for ii=3:length(headers)
    vr_data2p.(headers{ii})=interp1(vr_time,vr_data.(headers{ii}),time_vec_2p,'linear',0);
end
vr_data2p.timve_vec=time_vec_2p;
vr_data2p.valid_idx = valid_frames_2p;

%% downsample by number of layers
nLayers=2;
fn=fieldnames(vr_data2p);
for ii=1:length(fn)
    vr_data2p.(fn{ii})=vr_data2p.(fn{ii})(1:nLayers:end);
end



end

