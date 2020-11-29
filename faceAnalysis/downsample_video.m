function downsample_video(fn_in,fn_out)
vid = VideoReader(fn_in);
fn_ft = strrep(fn_in,'.avi','.ft');
ft = importdata(fn_ft);
vid_out = VideoWriter(fn_out);
ft_out = strrep(fn_out,'.avi','.ft');
vid_out.open()
nF=1;
oldFrame = [];

h = ones(3,3)/9;
frame_times_new=[];
try
while hasFrame(vid)
    video_1 = readFrame(vid);
    ft_1 = ft(nF,2);
    nF=nF+1;
    video_2 = readFrame(vid);
    ft_2 = ft(nF,2);
    nF=nF+1;
    fr_1 = imfilter(video_1,h);
    fr_2 = imfilter(video_2,h);
    new_frame = 0.5*fr_1+0.5*fr_2;
    
    vid_out.writeVideo(new_frame);
    frame_times_new(end+1)=0.5*ft_1+ft_2*0.5;
end

catch ME
    disp(ME.message)
    vid_out.close()
end
vid_out.close()
frame_vec = 0:numel(frame_times_new)-1;
frame_vec = [frame_vec',frame_times_new'];
writematrix(frame_vec,ft_out,'Delimiter','tab','FileType','text')


end
