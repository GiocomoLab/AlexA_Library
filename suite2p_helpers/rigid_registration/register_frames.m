function [shift_x,shift_y]=register_frames(data, template, perc_boundary)
%frame registration based on fft2
% register_frames(data, template, perc_boundary) will calculate the
% displacement values for x and y.
% 
% INPUT
%   data          - data stack
%   template      - reference image to which data stack is aligned (default is mean over stack)
%   perc_boundary - percentage of raws/columns cropped on image boundaries
%
% GK 01.05.2011
% ML 08.05.2014 docu

% default input values :
if nargin<2
    template=mean(data,3);
end

% user-defined values (hard-coded)


nbim_temp = 10; % template = mean over nbim_temp number of images
nbim_refresh = 50; %template refreshed 5 times per chunk

if nargin<3
    perc_boundary = 0.1; 
end


template=template-mean(template(:));
% determine how much of the images to use for alignment, the larger the
% boundary the less pixels are used for alignment and the faster the
% algorithm runs.
boundary=round(perc_boundary*max(size(data(:,:,1))));
if boundary >= size(data,1)
    boundary=round(perc_boundary*min(size(data(:,:,1))));
end
template=template(boundary+1:end-boundary,boundary+1:end-boundary);




% enabling parallel processing
%------------------------------
% if  matlabpool('size')==0
%     try
%         matlabpool open
%     end
% end
%------------------------------

high_pass_thresh=3;
low_pass_thresh=60;

upper_lim=round(prctile(template(:),95));
template(template>upper_lim)=upper_lim;
fft2_template = fft2(template); %fft2(mean(im2double(data_chunk(:,:,1:nbim_temp)),3));
fft2_template([1:high_pass_thresh size(fft2_template,1)-high_pass_thresh+2:size(fft2_template,1)],:)=0;
fft2_template(:,[1:high_pass_thresh size(fft2_template,2)-high_pass_thresh+2:size(fft2_template,2)])=0;
fft2_template(low_pass_thresh+2:end-low_pass_thresh,:)=0;
fft2_template(:,low_pass_thresh+2:end-low_pass_thresh)=0;
template_orig=template;



%disp(['Starting Chunk ' num2str(knd) ' of ' num2str(n_chunks) ])
%data_chunk=data(boundary+1:end-boundary,boundary+1:end-boundary, chunk_inds(knd):chunk_inds(knd+1)-1);
shifts=zeros(size(data,3),2);
for ind=1:size(data,3)
    if ~mod(ind,nbim_refresh)
        % refresh the template by aligning the last nbim_temp images
        temp_ali = shift_data(double(data(boundary+1:end-boundary,boundary+1:end-boundary,ind-nbim_temp:ind-1)),shifts(ind-nbim_temp:ind-1,1),shifts(ind-nbim_temp:ind-1,2));
        % ...and take the fft from the mean
        template=(double(mean(temp_ali,3))+template_orig)/2;
        upper_lim=round(prctile(template(:),95));
        template(template>upper_lim)=upper_lim;
        fft2_template = fft2(template);
        fft2_template([1:high_pass_thresh size(fft2_template,1)-high_pass_thresh+2:size(fft2_template,1)],:)=0;
        fft2_template(:,[1:high_pass_thresh size(fft2_template,2)-high_pass_thresh+2:size(fft2_template,2)])=0;
        fft2_template(low_pass_thresh+2:end-low_pass_thresh,:)=0;
        fft2_template(:,low_pass_thresh+2:end-low_pass_thresh)=0;
    end
    % calculates registration via DFT :
    [shifts(ind,1),shifts(ind,2)] = fftreg(fft2(double(data(boundary+1:end-boundary, boundary+1:end-boundary,ind))), fft2_template);
    
end


shift_x=shifts(:,1);
shift_y=shifts(:,2);























