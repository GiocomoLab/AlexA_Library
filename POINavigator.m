


flag =true;

while true

if flag
    
    POINav.hf=figure(1001);
    clf
    set(POINav.hf,'name','POINavigator','numberTitle','off');
    set(POINav.hf,'menubar','none');
    POINav.hf.UserData.FTA=5;

    POINav=POINav_createComponents(POINav);
    %mmfile.Format = {'int16' [1 16] 'header' ; 'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA'};
    %imageData = double(intmax('uint16')-mmfile.Data.chA);
     imageData = randn(256);
    
    POINav.hf.UserData.IMAGES=zeros([size(imageData),POINav.hf.UserData.FTA]);
    POINav.hf.UserData.COUNTER=0;
    flag = 0;

    colormap gray;      % use gray colormap
    
else
    %update graphics here
      %imageData = delta*imageData + (1-delta)*double(intmax('uint16')-mmfile.Data.chA);
        %ih.CData = mchA;
    newFrame=randn(256);
    POINav.hf.UserData.COUNTER=POINav.hf.UserData.COUNTER+1;
    idx=mod(POINav.hf.UserData.COUNTER,POINav.hf.UserData.FTA)+1;
    POINav.hf.UserData.IMAGES(:,:,idx)=newFrame;
    av_image=squeeze(mean(POINav.hf.UserData.IMAGES,3));
    if POINav.hf.UserData.OVERLAY
        %overlay with current template
        template_frame=POINav.Template_Ax_Image.CData;
        CData=imfuse(av_image,template_frame,'blend','Scaling','joint');
    else
        CData=av_image;
    end
        POINav.Live_Ax_Image.CData=CData;
        POINav.hf.UserData.currX=rand(1)*5-2.5;
        POINav.hf.UserData.currY=rand(1)*5-2.5;
        POINav.hf.UserData.currZ=rand(1)*5-2.5;

        
    
end

drawnow limitrate

end