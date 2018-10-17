


flag =true;

while true

if flag
    
    POINav.hf=figure(1001);
    clf
    set(POINav.hf,'name','POINavigator','numberTitle','off');
    set(POINav.hf,'menubar','none');
    POINav.hf.UserData.FTA=5;

    mmfile.Format = {'int16' [1 16] 'header' ; 'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA'};
    imageData = double(intmax('uint16')-mmfile.Data.chA);
    POINav=POINav_createComponents(POINav,size(imageData,1),size(imageData,2));

     %imageData = randn(256);
    
    POINav.hf.UserData.IMAGES=zeros([size(imageData),POINav.hf.UserData.FTA]);
    POINav.hf.UserData.COUNTER=0;
    flag = 0;

    colormap gray;      % use gray colormap
    
else
    %update graphics here
      %imageData = delta*imageData + (1-delta)*double(intmax('uint16')-mmfile.Data.chA);
        %ih.CData = mchA;
    newFrame=mmfile.Data.chA;
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
        %POINav.Template_Ax_Image.CData=squeeze(POINav.hf.UserData.templates(:,:,POINav.hf.UserData.selectedRow));
        POINav.hf.UserData.currX=mmfile.Data.header(10);
        POINav.hf.UserData.currY=mmfile.Data.header(11);
        POINav.hf.UserData.currZ=mmfile.Data.header(12);

        
    
end
mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
drawnow limitrate

end