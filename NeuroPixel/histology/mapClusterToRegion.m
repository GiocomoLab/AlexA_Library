% load excel db
%session_table = readtable('Z:\giocomo\attialex\NP_DATA\AA_session_summary.xlsx');
%histo_loc = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology';
%session_loc = 'Z:\giocomo\attialex\NP_DATA';
addpath(genpath('/home/users/attialex/spikes/'))
addpath(genpath('/home/users/attialex/AlexA_Library'))
session_table = readtable('/oak/stanford/groups/giocomo/attialex/NP_DATA/AA_session_summary.xlsx');
histo_loc = '/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/Histology';
session_loc = '/oak/stanford/groups/giocomo/attialex/NP_DATA';
%go through each row
nS=size(session_table,1);
for iS=11:nS
    animal = session_table.Mouse{iS};
    sessionID = session_table.SessionName{iS};
    recording_day =session_table.RecordingDay(iS);
%load sp file
    clear data
    data = load(fullfile(session_loc,sessionID));
    if ~session_table.Histology(iS)
        continue
    end
%locate and load probe table
    clear probe_table
    probe_table = load(fullfile(histo_loc,animal,'combined','processed',sprintf('probe_%d_region_table',recording_day)));
    if strcmp(sessionID,'AA1_190730_gaincontrast10_1')
        %because something is wrong with this one
        probe_table.probe_length = 3000;
    end

%assign a label to each cluster
    clear anatomy
    [cluster_region,cluster_parent,tip_distance,depth]=getClusterRegion(data.sp,probe_table.borders_table,probe_table.probe_length);
    anatomy.cluster_region = cluster_region;
    anatomy.cluster_parent = cluster_parent;
    anatomy.tip_distance = tip_distance';
    anatomy.depth = depth';
    save(fullfile(session_loc,sessionID),'anatomy','-append')
    set(gcf,'Name',sessionID)
    drawnow
end
    
    % go find rasterplot for each cluster

%put it in folder 'parent region'