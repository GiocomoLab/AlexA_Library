%% get session info

session_infos= dir(fullfile('Z:\giocomo\export\data\Projects\AlexA_NP\AA_200919_1*','session_info.csv'));


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
        if ~session_info.Histology(iR)
            disp('no histology')
            continue
        end
        
        if sagittal
            if hemisphere == 'R'
                extra='_right';
            elseif hemisphere =='L'
                extra = '_left';
            end
        else
            extra = '';
        end
        
        fn = fullfile(oak,session_info.Mouse{iR},'histology',[ session_info.SessionName{iR} '.csv']);
        if ~isfile(fn)
            fn=fullfile(oak,session_info.Mouse{iR},'histology',[ session_info.SessionName{iR} '_KS25.csv']);
        end
            anatomy = readtable(fn);
            
        glx_name=session_info.SessionName{iR};
        parts = strsplit(glx_name,'_');
        exp_day = parts{4};
        if startsWith(session_info.Mouse{iR},{'AA_200920_4','AA_200919_1'})
            disp('getting vr from excel files')
            file_table = readtable("Z:\giocomo\attialex\Histology\histology_october.xlsx");
            %rec_day = strcat(exp_day(5:6),'.',exp_day(3:4));
            rec_day = session_info.Date{iR}
            row_idx = startsWith(file_table.Mouse,session_info.Mouse{iR}) & startsWith(file_table.Date,rec_day);
            vr_files_exist = true;
            tmp_vr_files = file_table.SessionName(row_idx);
        else
            vr_files_exist = false;
        end
        

        %find vr sessions
        if ~vr_files_exist 
            vr_file = dir(fullfile(oak,session_info.Mouse{iR},['*' exp_day '*.log']));
        else
            vr_file = struct();
            for ivr=1:numel(tmp_vr_files)
                vr_file(ivr).folder = fullfile(oak,session_info.Mouse{iR});
                vr_file(ivr).name = strcat(tmp_vr_files{ivr},'.mat');
            end
        end
        
        for iVR=1:numel(vr_file)
            %find spmat file
            mat_name = strrep(vr_file(iVR).name,'.log','.mat');
            mat_file_path = fullfile(oak,session_info.Mouse{iR},'ks_data',glx_name,mat_name);
            if isfile(mat_file_path)
                save(mat_file_path,'anatomy','-append')
                copyfile(mat_file_path,'F:\Alex\new_2')
            else
                warning(['no .mat file for ' mat_name])
            end
            disp(isfile(mat_file_path))
        end
        
    end
end

% go find rasterplot for each cluster

%put it in folder 'parent region'
