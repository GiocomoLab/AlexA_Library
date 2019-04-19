function [pupilData] = analyzeFaceVideo(fn)
%fn='Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I1\videos\0416_baseline_1.avi';
vid = VideoReader(fn);
nFrames = vid.Duration*vid.FrameRate;
pupilData=nan(floor(nFrames)+100,3);

se = strel('disk',1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process relevant frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ii = 0;  % python is zero-indexed...this saves trouble later
while hasFrame(vid)
    ii=ii+1;
    frame = rgb2gray(readFrame(vid));
    %frame = rgb2gray(frame);
    % using the first VR frame, set the ROIs for pupil and whisk
    if ii == 1
        
        % make fig
        h = figure();
        ax = axes('Parent', h);
        imshow(frame, 'Parent', ax);
        title(ax, sprintf('Frame #%d', 1));
        set(gcf,'Units','normalized')
        ax.Position = [0,0,1,1];
        
        % draw pupil ROI, corners at eyelids
        disp('Ready to draw pupil ROI? (press any key to continue)')
        pause
        [~, eye_corners] = makeROI(frame, h);
        disp('Ready to identify pupil? (press any key to continue)')
        pause
        ex = frame(eye_corners(1):eye_corners(2),eye_corners(3):eye_corners(4));
        imshow(ex,'Parent',ax)
        [~,center_corners] = makeROI(ex, h);
        level=thresh_tool(ex);
        clf
        [bwm,xi,yi]=roipoly(ex');
        disp('Processing video data (this may take awhile)...')
    else
        ex=frame(eye_corners(1):eye_corners(2),eye_corners(3):eye_corners(4));
    end
    
    %threshold image
    bw=ex<level;
    bw=bw&bwm';
    bw1=imdilate(bw,se);
    L=bwlabel(bw1);
    pup_center = L(center_corners(1):center_corners(2),center_corners(3):center_corners(4));
    id=max(pup_center(:));
    bw2=bwconvhull(L==id);
    outline=bwboundaries(bw2,'noholes');
    if id>0
        idx = inpolygon(outline{1}(:,1),outline{1}(:,2),xi,yi);
        a=CircleFitByPratt(outline{1}(idx,:));
        th = 0:pi/50:2*pi;
        xunit = a(3)* cos(th) + a(1);
        yunit = a(3) * sin(th) + a(2);
        pupilData(ii,:)=a;
    end
    if mod(ii,200)==1
        cla
        subplot(2,1,1)
        imagesc(bw2');
        subplot(2,1,2)
        imagesc(ex')
        hold on
        plot(xunit,yunit)
        plot(outline{1}(idx,1),outline{1}(idx,2),'r')
        %pause
    end
end
pupilData=pupilData(1:ii,:);
try
save([fn(1:end-4) '_pupilData.mat'],'pupilData');
catch ME
    disp('something went wrong while saving')
end
end
