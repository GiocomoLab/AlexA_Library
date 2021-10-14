function [cluster_table]=getClusterRegion_v2(tbl_row,borders_table,trajectory,metrics)

if ~isfield(metrics,'tip_distance')
ycoords = [];
for ii=1:192
    tmp = (ii-1)*20;
    ycoords = cat(1,ycoords,[tmp;tmp]);
end
tip_distance=ycoords(metrics.peak_channel+1);
else
tip_distance = metrics.tip_distance;
end
vars={'x','y','z'};
for iVar = 1:3
    var = vars{iVar};
    dd=max(tip_distance)-min(tip_distance);
    tmp = (trajectory.(var)(1)-trajectory.(var)(2))/dd * tip_distance*10 + trajectory.(var)(2)*10; %because first pos is brain exit
    pos.(var)=tmp;
    %pos.(var)=diff(trajectory.(var))/dd*tip_distance*10+trajectory.(var)(1)*10; %from pixel to mm
end


%%
cluster_position = tbl_row.ProbeDepth*1000-tip_distance;
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
if ~isrow(cluster_region)
cluster_region = cluster_region';
end
if ~isrow(cluster_parent)
cluster_parent = cluster_parent';
end
cluster_table = table(metrics.cluster_id',cluster_region',cluster_parent',pos.x',pos.y',pos.z',depth','VariableNames',{'cluster_id','cluster_region','cluster_parent','xpos','ypos','zpos','depth'});
end