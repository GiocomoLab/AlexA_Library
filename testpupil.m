cx=40;
cy=90;
figure
hold on
th1 =20;
th2 = 90;
se = strel('disk',1);
ex=idata(400:520,440:580,:);
for ii=1:1000
bw=ex(:,:,ii)<th1;
bw1=imdilate(bw,se);
%[outline,L]=bwboundaries(bw1,'noholes');
L=bwlabel(bw1);
id=L(cx,cy);
subplot(2,1,1)
bw3=bwconvhull(L==id);
imagesc(bwconvhull(L==id)');
[outline,L]=bwboundaries(bw3,'noholes');
%imagesc(L');
if id>0
a=CircleFitByPratt(outline{1});
th = 0:pi/50:2*pi;
xunit = a(3)* cos(th) + a(1);
yunit = a(3) * sin(th) + a(2);


subplot(2,1,2)
title(sprintf('%d',ii))

imagesc(ex(:,:,ii)')
hold on
plot(xunit, yunit);

end
pause
cla
end