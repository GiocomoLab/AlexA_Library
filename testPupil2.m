mov=VideoReader('Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels\I1\videos\0416_baseline_1.avi');
nFrames = mov.Duration*mov.FrameRate;
cx=40;
cy=90;
figure
hold on
th1 =20;
th2 = 90;
se = strel('disk',1);
ii=0;
pupilData=nan(floor(nFrames)+100,3);
while hasFrame(mov)
    %read frame and extract roi
    fr = rgb2gray(readFrame(mov));
    ii=ii+1;
    ex=fr(400:520,440:580);
    %threshold image
    bw=ex<th1;
    bw1=imdilate(bw,se);
    L=bwlabel(bw1);
    id=L(cx,cy);
    subplot(2,1,1)
    bw2=bwconvhull(L==id);
    [outline,L]=bwboundaries(bw2,'noholes');
    if id>0
        a=CircleFitByPratt(outline{1});
        th = 0:pi/50:2*pi;
        xunit = a(3)* cos(th) + a(1);
        yunit = a(3) * sin(th) + a(2); 
        pupilData(ii,:)=a;
    end
    if mod(ii,200)==0
        cla
        subplot(2,1,1)
        imagesc(bw2');
        subplot(2,1,2)
        imagesc(ex')
        hold on
        plot(xunit,yunit)
        
        
    
    end
end