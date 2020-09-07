startidx=1:30:591;
MM=zeros(512,796,length(startidx));
for ii=1:length(startidx)
    starti=startidx(ii);
    stopidx=min([starti+30,591]);
    MM(:,:,ii)=mean(XX(:,:,starti:stopidx),3);
    imwrite(uint8(mat2gray(squeeze(MM(:,:,ii)))*255),'test2.tif','WriteMode','append')
end