function  showXcorrResults(location,region,gains)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
figure('Position',[405         608        1355         370])
cm=lines(numel(gains));
leg={};
for ii=1:numel(gains)
    try
    tmp = fullfile(location,sprintf('allData_%s_%s_100.mat',region,gains{ii}));
    data = load(tmp);
    
subplot(1,2,1)
boundedline(data.output.X-3200,nanmean(data.output.Y),nanstd(data.output.Y)/sqrt(size(data.output.Y,1)),'alpha','cmap',cm(ii,:))
hold on
axis tight
ylim([0 0.8])
leg{end+1}=gains{ii};
subplot(1,2,2)

boundedline(data.output.X-3200,nanmean(data.output.S*2),nanstd(data.output.S*2)/sqrt(size(data.output.Y,1)),'alpha','cmap',cm(ii,:))
    axis tight
    ylim([-2 6])
    catch ME
    end
end
subplot(1,2,1)
legend(leg)
title(region)
ylabel('Peak XCorr')
xlabel('Distance from Gain onset')
subplot(1,2,2)
legend(leg)
title(region)
ylabel('Shift [cm]')
xlabel('Distance from Gain onset')
end

