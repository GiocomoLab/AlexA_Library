dat08 = load('Z:\giocomo\attialex\Images\xcorrv5\allData_MEC_0.8_100.mat');
dat05 = load('Z:\giocomo\attialex\Images\xcorrv5\allData_MEC_0.5_100.mat');
figure
subplot(1,2,1)
boundedline(dat08.output.X-4000,nanmean(dat08.output.Y),nanstd(dat08.output.Y)/sqrt(size(dat08.output.Y,1)),'alpha')
hold on
boundedline(dat05.output.X-4000,nanmean(dat05.output.Y),nanstd(dat05.output.Y)/sqrt(size(dat05.output.Y,1)),'r','alpha')
title('100% contrast')
legend({'0.8','0.5'})
subplot(1,2,2)
boundedline(dat08.output.X-4000,nanmean(dat08.output.S),nanstd(dat08.output.S)/sqrt(size(dat08.output.Y,1)),'alpha')
hold on
boundedline(dat05.output.X-4000,nanmean(dat05.output.S),nanstd(dat05.output.S)/sqrt(size(dat05.output.Y,1)),'r','alpha')

%%


dat08 = load('Z:\giocomo\attialex\Images\xcorrv5\allData_MEC_0.8_10.mat');
dat05 = load('Z:\giocomo\attialex\Images\xcorrv5\allData_MEC_0.5_10.mat');
figure
subplot(1,2,1)
boundedline(dat08.output.X,nanmean(dat08.output.Y),nanstd(dat08.output.Y)/sqrt(size(dat08.output.Y,1)),'alpha')
hold on
boundedline(dat05.output.X,nanmean(dat05.output.Y),nanstd(dat05.output.Y)/sqrt(size(dat05.output.Y,1)),'r','alpha')
title('10% contrast')

legend({'0.8','0.5'})
subplot(1,2,2)
boundedline(dat08.output.X-4000,nanmean(dat08.output.S),nanstd(dat08.output.S)/sqrt(size(dat08.output.Y,1)),'alpha')
hold on
boundedline(dat05.output.X-4000,nanmean(dat05.output.S),nanstd(dat05.output.S)/sqrt(size(dat05.output.Y,1)),'r','alpha')
