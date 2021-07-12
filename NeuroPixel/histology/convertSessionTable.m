session_table = readtable('Z:\giocomo\attialex\Histology\session_february.xlsx');
root = 'Z:\giocomo\export\data\Projects\AlexA_NP';
%%
old_mouse = '';
old_data = '';
for iS = 1:size(session_table,1)
    current_mouse = session_table.Mouse{iS};
    if ~isempty(session_table.Mouse{iS})
        if ~strcmp(current_mouse,old_mouse)
            %fprintf('new mouse: %s \n',current_mouse)
            old_mouse = current_mouse;
            Mouse = {};
            Date = {};
            Hemisphere = [];
            ProbeDepth = [];
            SessionName= {};
            Histology=[];
            Sagittal = [];
            RecordingDay =[];
            idx = startsWith(session_table.Mouse,current_mouse);
            dates = unique(session_table.Date(idx));
            for iDate = 1:numel(dates)
                current_date = dates{iDate};
                combined_idx = find(startsWith(session_table.Mouse,current_mouse) & startsWith(session_table.Date,current_date),1);
                [dd]=strsplit(session_table.Date{combined_idx},'.');
                folder = dir(fullfile(root,session_table.Mouse{combined_idx},[session_table.Mouse{combined_idx}, '*', dd{3}(1:2), dd{2},dd{1},'*_g0']));

                if numel(folder)==1
                    fo = folder(1).name;
                else
                    folder = dir(fullfile(root,session_table.Mouse{combined_idx},'neuropixels_data',[session_table.Mouse{combined_idx}, '*', dd{3}(1:2), dd{2},dd{1},'*_g0']));
                    if numel(folder)==1
                        fo=folder(1).name;
                    else
                    fo = ' ';
                    end
                end
                Mouse = cat(1,Mouse,session_table.Mouse{combined_idx});
                Date = cat(1,Date,session_table.Date{combined_idx});
                Hemisphere = cat(1,Hemisphere,session_table.Hemisphere{combined_idx});
                ProbeDepth = cat(1,ProbeDepth,session_table.ProbeDepth(combined_idx));
                Histology = cat(1,Histology, session_table.Histology(combined_idx));
                SessionName = cat(1,SessionName,fo);
                Sagittal = cat(1,Sagittal,session_table.Sagittal(combined_idx));
                RecordingDay = cat(1,RecordingDay,session_table.RecordingDay(combined_idx));
            end
            tbl = table(Mouse,Date,Hemisphere,Histology,ProbeDepth,SessionName,Sagittal,RecordingDay);
            fn = fullfile(root,Mouse{1},'session_info.csv');
            writetable(tbl,fn);
        end
    end
    
end
    
    