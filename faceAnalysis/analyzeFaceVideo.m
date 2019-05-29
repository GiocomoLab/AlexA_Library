function [pupilData] = analyzeFaceVideo(fn)
%fn='Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I1\videos\0416_baseline_1.avi';
vid = VideoReader(fn);
nFrames = vid.Duration*vid.FrameRate;
pupilData=nan(floor(nFrames)+100,3);

se = strel('disk',1);

if exist([fn(1:end-4) '_pupilData.mat'],'file')
    fprintf('now working on %s',fn)
    s=input('PupilData file for this already exists, skip [0/1]:');
    if s==1
        return
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process relevant frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vid.CurrentTime = 10;

ii = 0;  % python is zero-indexed...this saves trouble later
h2=figure('Name',fn);
ax2=axes('Parent',h2);
while hasFrame(vid)
    ii=ii+1;
    frame = rgb2gray(readFrame(vid));
    %frame = rgb2gray(frame);
    % using the first VR frame, set the ROIs for pupil and whisk
    if ii == 1
        ct=vid.CurrentTime;
        vid.CurrentTime = 10;
        template = rgb2gray(readFrame(vid));
        vid.CurrentTime = ct;
        
        % make fig
        h = figure('Name',fn);
        ax = axes('Parent', h);
        imshow(template, 'Parent', ax);
        title(ax, sprintf('Frame #%d', 1));
        set(gcf,'Units','normalized')
        ax.Position = [0,0,1,1];
        
        [~, eye_corners] = makeROI(frame, h);
        disp('Ready to identify pupil? (press any key to continue) (subselect valid pupil area) \n only point within area will be considered for circle fitting')
        pause
        ex = template(eye_corners(1):eye_corners(2),eye_corners(3):eye_corners(4));
        imshow(ex,'Parent',ax)
        [~,center_corners] = makeROI(ex, h);
        level=thresh_tool(ex);
        clf
        [bwm,xi,yi]=roipoly(ex');
        close(h);
        disp('Processing video data (this may take awhile)...')
        figure(h2)
        axes(ax2)
%     else
%         ex=frame(eye_corners(1):eye_corners(2),eye_corners(3):eye_corners(4));
    end
    ex=frame(eye_corners(1):eye_corners(2),eye_corners(3):eye_corners(4));
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
    if mod(ii,100)==1
        cla
        subplot(2,1,1)
        %imagesc(bw2');
        plot(pupilData(1:ii,1))
        subplot(2,1,2)
        imagesc(ex')
        hold on
        plot(xunit,yunit)
        plot(outline{1}(idx,1),outline{1}(idx,2),'r')
        title(sprintf('%d',ii))
        
         %keyboard
    end
end
pupilData=pupilData(1:ii,:);
try
    figure('Name',fn)
    plot(pupilData(:,3))
    title(sprintf('using level %.1f',level))
save([fn(1:end-4) '_pupilData.mat'],'pupilData','level');
catch ME
    disp('something went wrong while saving')
end
end
