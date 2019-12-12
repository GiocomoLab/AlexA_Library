params = readtable('UniversalParams.xlsx');

figure

for cellIDX = 1:size(spatialMap,1)
    
    if ~startsWith(reg(cellIDX),'VISp')
        continue
end
    
    
    
    SMS=squeeze(spatialMap(cellIDX,:,:));
    mS=mean(SMS,2);
    [ma,mi]=max(mS);
    
    posWindow=[mi*2-30 mi*2-27];
    data.posWindow = posWindow;

    
    subplot(1,4,1)
    plot(mS)
    hold on
    plot(mi,ma,'ro')
    subplot(1,3,2)
%     imagesc(SMS')
%     ylim([1 16])
%     xlim([mi-10 mi+1])

    hold on
    trialsorted = plotRasterSpeedSort(data,params,[],good_cells(cellIDX),trials,subplot(1,4,2));
        xlim([mi*2-30 mi*2+1])

    subplot(1,4,3)
     spike_id=data.sp.clu==good_cells(cellIDX);
     spike_t = data.sp.st(spike_id);
     [~,~,spike_idx] = histcounts(spike_t,data.post);
     scatter(data.posx(spike_idx),data.trial(spike_idx)-4,2)
    %ylim([min(trials)-1 max(trials+1)])
    xlim([mi*2-30 mi*2+1])
    ylim([1 16])
    %xlim([mi-10 mi+1])
    subplot(1,4,4)
    tr = trialsorted(trials);
    mS=mean(SMS(:,tr<numel(trials)/4),2);
    
    plot(mS)
    hold on
    mS=mean(SMS(:,tr>numel(trials)/4*3),2);
    plot(mS)
    xlim([mi-30 mi+1])

    pause
    clf
end