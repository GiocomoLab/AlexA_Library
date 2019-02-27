function [dx,dy,template,registered_data] = rigid_registration(data)
template_frames=11:30;
fprintf('...registering...\n')
[dx,dy]=register_frames(data,mean(data(:,:,template_frames),3));
fprintf('...shifting data and correcting line shift...\n')
registered_data=shift_data(data,dx,dy);
registered_data=correct_line_shift(registered_data,mean(registered_data,3));
template=mean(registered_data,3);

end
