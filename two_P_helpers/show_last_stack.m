function varargout=show_last_stack(file_name)
% show_last_stack loads .bin file and displays raw 2P imaging data in a
% movie.
%
% Loads the last stack found in imaging temp directory or stack specified by ID
%
% OPTIONAL inputs:
% ch_combined (default=1) - combine channels for two channel data (needs stackID)
% splitdata (default=1)   - splitdata if multilayer stack
% calcMM (default=0)      - plot mismatch response based on full-frame
%                           averaged signal. Only works with single channel
%                           and splitdata
%e.g.: show_last_stack();
%e.g.: show_last_stack(12345);
%
%Control for movie as for tiff_stack
% 'shift+arrow up/down' - adjust brightness
% 'arrow up/down' - adjust movie speed
% 'space tab' - pause movie
% 't' - create average image of stack

%

if nargin<1
    file_name=find_file_name();
end


%load info file m file
load([file_name '.mat']);
%get number of frames
stop=max(info.frame);

data=sbxread(file_name,0,stop);

view_stack(squeeze(data(1,:,:,:)));

if nargout>0 %prevent console flooding 
    varargout{1}=data;
    varargout{2}=data1;
end

end

function check_zdrift(data)
[dx,dy] = register_frames(data,mean(data(:,:,20:30),3));
data = shift_data(data,dx,dy);
im1 = mean(data(:,:,1:floor(end/2)),3);
im2 = mean(data(:,:,ceil(end/2+1):end),3);
im3(:,:,1) = ntzo(im1);
im3(:,:,2) = ntzo(im2);
im3(:,:,3) = zeros(size(im1));
figure;
colormap gray
subplot(1,3,1); imagesc(imadjust(ntzo(im1),[0 0.4]));axis off;title('first 50 frames')
subplot(1,3,2); imagesc(imadjust(ntzo(im2),[0 0.4]));axis off;title('last 50 frames')
subplot(1,3,3); imagesc(im3);axis off;title('combined')
end

function fn=find_file_name()
aa=importdata('D:\\datalog.txt');

fn=strtrim(aa{end});


end