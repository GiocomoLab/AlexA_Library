% load excel db
%session_table = readtable('Z:\giocomo\attialex\NP_DATA\AA_session_summary.xlsx');
%histo_loc = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology';
%session_loc = 'Z:\giocomo\attialex\NP_DATA';
addpath(genpath('C:\code\spikes'))
addpath(genpath('/home/users/attialex/AlexA_Library'))
session_table = readtable('Z:\giocomo\attialex\Histology\histology_october.xlsx');
%histo_loc = '/oak/stanford/groups/giocomo/export/data/Projects/AlexA_NP/Histology';
session_loc = 'F:\Alex\matfiles_new';
%go through each row
nS=size(session_table,1);
oak = 'Z:\giocomo\export\data\Projects\AlexA_NP';
%%
for iS=51:nS
    if isempty(session_table.Mouse{iS})
        continue
    end
     if ~session_table.Histology(iS)
        continue
    end
    animal = session_table.Mouse{iS};
    sessionID = session_table.SessionName{iS};
    recording_day =session_table.RecordingDay(iS);
%load sp file
    clear data
    data = load(fullfile(session_loc,sessionID));
   
%locate and load probe table
    histo_loc = fullfile(oak,animal,'histology');
    clear probe_table
    hemisphere = session_table.Hemisphere{iS};
    sagittal = session_table.Sagittal(1);
    
    if sagittal
        if hemisphere == 'R'
            extra='_right';
        elseif hemisphere =='L'
            extra = '_left';
        end
    else
        extra = '';
    end
    probe_table = load(fullfile(histo_loc,['combined',extra],['processed'],sprintf('probe_%d_region_table.mat',recording_day)));
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
    anatomy.probe_trajectory = probe_table.trajectory;
    save(fullfile(session_loc,sessionID),'anatomy','-append')
    new_loc = fullfile(oak,animal,'ks_data');
    copyfile(fullfile(session_loc,[sessionID '.mat']),new_loc)
    set(gcf,'Name',sessionID)
    saveas(gcf,fullfile(new_loc,[sessionID '_histo.png']))
    drawnow
end
    
    % go find rasterplot for each cluster

%put it in folder 'parent region'