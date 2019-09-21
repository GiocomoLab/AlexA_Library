%session_table = readtable('Z:\giocomo\attialex\NP_DATA\AA_session_summary.xlsx');
%histo_loc = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology';
%session_loc = 'Z:\giocomo\attialex\NP_DATA';
%%
session_table = readtable('/oak/stanford/groups/giocomo/attialex/NP_DATA/AA_session_summary.xlsx');
histo_loc = '/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/Histology';
session_loc = '/oak/stanford/groups/giocomo/attialex/NP_DATA';
%go through each row
nS=size(session_table,1);
%%
for iS=6:nS
    animal = session_table.Mouse{iS};
    sessionID = session_table.SessionName{iS};
    recording_day =session_table.RecordingDay(iS);
    recording_depth = session_table.ProbeDepth(iS);
%load sp file
    clear data
    clear anatomy
    data = load(fullfile(session_loc,sessionID),'anatomy');
    if ~session_table.Histology(iS)
        continue
    end
    figure
    subplot(1,2,1)
    plot(data.anatomy.tip_distance,'.')
    yline(recording_depth*1000)
    
    d1=recording_depth*1000-max(data.anatomy.tip_distance);
    d2=min(data.anatomy.depth);
    title(sessionID)
    xlabel(sprintf('D1: %.2f, D2: %.2f',d1,d2))
    cluster_position = data.anatomy.depth;
    probe_table = load(fullfile(histo_loc,animal,'combined','processed',sprintf('probe_%d_region_table',recording_day)));
    borders_table = probe_table.borders_table;
    subplot(1,2,2)
    x = randn(size(cluster_position));
    plot(x,cluster_position,'.')
    if d2>200
        cp_shifted = cluster_position-d2+200;
        hold on
        plot(x,cp_shifted,'.')
    end
    midY=(borders_table.upperBorder+borders_table.lowerBorder)/2;
    acr = borders_table.acronym;
    set(gca, 'YTick', midY, 'YTickLabel', acr);
    set(gca, 'YDir','reverse');
    
    if d2>200
       cp_shifted = cluster_position-d2+200;
       lower= borders_table.lowerBorder;
    upper = borders_table.upperBorder;
    acr={};
    parent = {};
    for iC =1:length(cp_shifted)
        idx = find(cp_shifted(iC)>=upper & cp_shifted(iC)<lower);
        if ~isempty(idx)
        acr{iC}=borders_table.acronym{idx};
        parent{iC}=borders_table.parent_name{idx};
        else
            acr{iC}='';
            parent{iC}='';
        end
    end

anatomy = data.anatomy;
anatomy.depth_shifted = cp_shifted;
anatomy.region_shifted = acr;
anatomy.parent_shifted = parent;
        save(fullfile(session_loc,sessionID),'anatomy','-append')

    end
    

    drawnow
    
end
%%