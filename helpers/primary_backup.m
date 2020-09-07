function primary_backup(mouse_name,exp_name,varargin)
% copies files from local machine to server
% usage: primary_backup('F2','1015_contrasttrack_gainchanges_2')
% optional parameters:  'root_dir' location of 'mouse_name folder' to backup , e.g. F:
%                       'server_root' of 'mouse_name folder on backup
%                       location, e.g.
%                       Z:\Projects\ContrastExperiment_neuropixels
%                       'relevant_endings': cell array of filetypes to copy
%                       TODO: copy all files if empty cell array is passed
% expects root\mouse\expName\, will create a folder on backup dir:
% serverroot\mouse\expName (if not already existign) and copies specified
% filetypes to that location


p = inputParser;
   default_root_dir = 'F:';
   default_server_root = 'Z:\Projects\ContrastExperiment_neuropixels';
   default_endings= {'*.npy','*.mat','*.m','*.py','*.tsv'};
   
   addParameter(p,'root_dir',default_root_dir);
   addParameter(p,'server_root',default_server_root);
   addParameter(p,'relevant_endings',default_endings);
   
   parse(p,varargin{:});

relevant_endings=p.Results.relevant_endings;
root_dir = p.Results.root_dir;
server_root = p.Results.server_root;
if ~isdir(fullfile(root_dir,mouse_name,exp_name))
    disp(fullfile(root_dir,mouse_name,exp_name))
    error('mouse and/or exp not found on specified directory')
end

if ~isdir(fullfile(server_root,mouse_name))
    disp(fullfile(server_root,mouse_name))
    error('mouse not found on server')
end

backup_files={};
size=0;
for iE = 1:length(relevant_endings)
    fn=fullfile(root_dir,mouse_name,exp_name,relevant_endings{iE});
    ff = dir(fn);
    
    
    for iF=1:length(ff)
        backup_files{end+1}=fullfile(ff(iF).folder,ff(iF).name);
        size=size+ff(iF).bytes;
    end
end
fprintf('about to copy %d files, %.3f MB \n',length(backup_files),size/1000000)
%%
backup_copy(backup_files,fullfile(server_root,mouse_name,exp_name))

%%

fnames_local = dir(fullfile(root_dir,mouse_name,exp_name));

%// Extract the other .tgz file to the same location
fnames_server = dir(fullfile(server_root,mouse_name,exp_name));

%// Use setdiff to compare the files that were in one but not the other
in_local_but_not_server = setdiff({fnames_local.name}, {fnames_server.name});
fprintf('Not on server:\n')
fprintf(1, '%s \n ', in_local_but_not_server{:})
