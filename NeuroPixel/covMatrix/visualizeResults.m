showXcorrResults('Z:\giocomo\attialex\Images\xcorrv7','MEC',{'1.0','0.8','0.5'})
    
%%
shiftfig = figure();
xtxfig = figure();
levels = {'50','20','10'};
cm=lines(numel(levels));
region = 'MEC';
for iL=1:numel(levels)
    figure(shiftfig)
data =load(['Z:\giocomo\attialex\Images\xcorrv_CONTRAST2\allData_' region '__' levels{iL} '.mat']);
subplot(1,2,1)
boundedline(data.output.X-2400,nanmean(data.output.Y),nanstd(data.output.Y)/sqrt(size(data.output.Y,1)),'alpha','cmap',cm(iL,:))
title('XCorr peak values')
ylabel('XCorr')
xlabel('Distance from contrast change [cm]')
subplot(1,2,2)
boundedline(data.output.X-2400,nanmean(data.output.S*2),nanstd(data.output.S*2)/sqrt(size(data.output.Y,1)),'alpha','cmap',cm(iL,:))
title('XCorr Shifts')
ylabel('Shift [cm]')
xlabel('Distance from contrast change [cm]')
figure(xtxfig)
subplot(2,numel(levels),iL)
imagesc(data.output.XTX,[0 .8])
set(gca,'XTickLabel',[],'YTickLabel',[])
xline(1200,'r')
yline(1200,'r')
axis image
title(levels{iL})
colorbar
subplot(2,numel(levels),numel(levels)+iL)
imagesc(nanmean(data.output.YYT,3),[0. .5])
set(gca,'XTickLabel',[],'YTickLabel',[])
colorbar
title(levels{iL})

axis image
end
figure(shiftfig)
subplot(1,2,1)
legend(levels)
%%
figure('Position',[680   558   926   420])
plot(output.Y(7:10,:)')
set(gca,'XTick',[5:10:200])
set(gca,'XTickLabel',[-6:1:-1 1:4 1:6])
xline(6*10+0.5,'r')
xline(10*10+0.5,'r')
ylabel('Peak Xcorr')
xlabel('Trial #')
