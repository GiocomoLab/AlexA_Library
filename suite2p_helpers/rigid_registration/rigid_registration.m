function [dx,dy,template,tmpdat] = rigid_registration(data)
template_frames=11:30;
fprintf('...registering...\n')
[dx,dy]=register_frames(data,mean(data(:,:,template_frames),3));
fprintf('...shifting data and correcting line shift...\n')
tmpdat=shift_data(data,dx,dy);
tmpdat=correct_line_shift(tmpdat,mean(tmpdat,3));
template=mean(tmpdat,3);

end
