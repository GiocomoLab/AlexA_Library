function [roi, roi_corners] = makeROI(image, fig_handle, varargin)

% IL 7/20/18
% draw an ROI on a figure
%
% input image matrix, figure handle, and color of ROI
% outputs the ROI as a normalized, grayscale image matrix
% also outputs row and column indices to index into image matrix and
% extract ROI
% ex. roi = image(roi_corners(1):roi_corners(2),roi_corners(3):roi_corners(4));

% draw ROI
if nargin < 3
    figure(fig_handle)
    k = waitforbuttonpress;
    rect_pos = rbbox;
    annotation('rectangle',rect_pos,'Color','red')
    x0 = rect_pos(1) * size(image,2);
    x1 = rect_pos(3) * size(image,2);    
    x0 = round(x0); x1 = round(x1);
    y0 = size(image,1) - (rect_pos(2) * size(image,1));
    y1 = rect_pos(4) * size(image,1);
    y0 = round(y0); y1 = round(y1);
    varargin = [(y0-y1) y0 x0 (x0+x1)];    
end

if nargin == 3
    roi_corners = varargin{1};
else
    roi_corners = varargin;
end

% index into image matrix
roi = image(roi_corners(1):roi_corners(2),roi_corners(3):roi_corners(4));
roi = im2double(roi);
roi = mat2gray(roi,[min(min(roi)) max(max(roi))]); % normalize
    
end