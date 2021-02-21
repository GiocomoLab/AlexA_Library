function [cluster_region,cluster_parent,tip_distance,depth]=getClusterRegion(sp,borders_table,probe_length)

if isfield(sp,'waveform_metrics')
    ycoords=[];

for ii=1:192
    tmp = (ii-1)*20;
    ycoords = cat(1,ycoords,[tmp;tmp]);
end
    tip_distance=ycoords(sp.waveform_metrics.peak_channel+1);
else
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
        templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
    %%
    tip_distance=zeros(length(sp.cgs),1);
    for iC=1:length(sp.cgs)
        tip_distance(iC)=nanmean(spikeDepths(sp.clu==sp.cids(iC)));
    end
end
    %%
    cluster_position = probe_length-tip_distance;
    depth=cluster_position;
    if nnz(depth<0)>0
        warning('spikes outside brain')
%         %keyboard
%         aa=find(sp.cgs==2 & depth'<0); %because we want to avoid noise or mua clusters for remapping this
%         if isempty(aa)
%             offset = max(tip_distance);
%         else
%             offset = max(tip_distance(aa));
%         end
%             td_norm = (tip_distance-min(tip_distance))/(offset-min(tip_distance));
%         td_norm = td_norm*probe_length;
%         cluster_position = probe_length-td_norm;
%         cluster_position(cluster_position<0)=0;
%         depth = cluster_position;
    end
    
    %% plot clusters as a function of depth, with labels
    figure('Position',[680   208   389   770]);
    plot(randn(size(cluster_position)),cluster_position,'.')
    midY=(borders_table.upperBorder+borders_table.lowerBorder)/2;
    acr = borders_table.acronym;
    set(gca, 'YTick', midY, 'YTickLabel', acr);
    set(gca, 'YDir','reverse');
    yyaxis('right')
    set(gca,'YTick',borders_table.upperBorder)
    set(gca,'YDir','reverse')
    
    set(gca,'YLim',[0, max(borders_table.lowerBorder)])
    yyaxis('left')
    set(gca,'YLim',[0, max(borders_table.lowerBorder)])

    %% region labels for each cluster
    lower= borders_table.lowerBorder;
    upper = borders_table.upperBorder;
    acr={};
    parent = {};
    for iC =1:length(tip_distance)
        idx = find(cluster_position(iC)>=upper & cluster_position(iC)<lower);
        if ~isempty(idx)
        acr{iC}=borders_table.acronym{idx};
        parent{iC}=borders_table.parent_name{idx};
        else
            acr{iC}='';
            parent{iC}='';
        end
    end
cluster_region = acr;
cluster_parent = parent;
    
end