function backup_copy(files_to_copy,targetDir)



if ~isdir(targetDir)
    mkdir(targetDir)
end
no_backup_files = length(files_to_copy);

backup_files = {};

for ind = 1:no_backup_files
    [folder, fn,ending]=fileparts(files_to_copy{ind});
    backup_files{end+1} = fullfile(targetDir,[fn ending]);
end

warnings={};
for ind = 1:no_backup_files
    % check whether directory exist
    fprintf('copying %d of %d \n',ind,no_backup_files)
    if ~exist(backup_files{ind},'file')
%         disp(['-> ' sources_files{ceil(ind/no_backup_loc)} ' >> ' backup_files{ind}]);
        [success(ind),mtext] = copyfile(files_to_copy{ind},backup_files{ind});
    else
        warnings{end+1} = ([backup_files{ind} ' - FILE ALREADY EXISTS - skipping.']);
        warning(warnings{end});
        success(ind) = 1;
    end
    if ~success(ind)
        warnings{end+1} = strtrim([backup_files{ind} ' - Some problem with copying: ' mtext]);
        warning(warnings{end});
        success(ind) = 0;
    end
end