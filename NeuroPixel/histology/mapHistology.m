%% get session info

session_infos= dir(fullfile('Z:\giocomo\export\data\Projects\AlexA_NP\AA_210114*','session_info.csv'));


oak = 'Z:\giocomo\export\data\Projects\AlexA_NP';
%%  
for iS=1:size(session_infos,1)
    session_info = readtable(fullfile(session_infos(iS).folder,session_infos(iS).name),'Delimiter',',');
    for iR = 1:size(session_info,1)
        %locate and load probe table
        histo_loc = fullfile(oak,session_info.Mouse{iR},'histology');
        clear probe_table
        hemisphere = session_info.Hemisphere{iR};
        sagittal = session_info.Sagittal(iR);
        
        if sagittal
            if hemisphere == 'R'
                extra='_right';
            elseif hemisphere =='L'
                extra = '_left';
            end
        else
            extra = '';
        end
        probe_table = load(fullfile(histo_loc,['combined',extra],['processed'],sprintf('probe_%d_region_table.mat',session_info.RecordingDay(iR))));
        metrics_file = fullfile(oak,session_info.Mouse{iR},'ks_data',session_info.SessionName{iR},'waveform_metrics.csv');
        %metrics_file = fullfile('H:\catGT2',['catgt_' session_info.SessionName{iR}],[session_info.SessionName{iR} '_imec0'],'imec0_ks2','waveform_metrics.csv');
        
        metrics = readtable(metrics_file);
        cluster_table=getClusterRegion_v2(session_info(iR,:),probe_table.borders_table,probe_table.trajectory,metrics);
        
        set(gcf,'Name',session_info.SessionName{iR})
        sn_img = fullfile(oak,session_info.Mouse{iR},'histology',[ session_info.SessionName{iR} '_KS25_histo.png']);
        fn = fullfile(oak,session_info.Mouse{iR},'histology',[ session_info.SessionName{iR} '_KS25.csv']);
        saveas(gcf,sn_img)
        drawnow
        writetable(cluster_table,fn);
    end
end

% go find rasterplot for each cluster

%put it in folder 'parent region'
