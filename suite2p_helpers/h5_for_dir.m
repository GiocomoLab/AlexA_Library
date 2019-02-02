function h5_for_dir(input_dir)

sbx_files = dir(fullfile(input_dir,'*.sbx'));

for ii=1:length(sbx_files)
    fprintf('Now working on: %s \n',sbx_files(ii).name)
    [~,fn,~]=fileparts(sbx_files(ii).name);
    if ~exist([fn, '.h5'])
        sbx2h5(fullfile(input_dir,sbx_files(ii).name))
        fprintf('done /n')
        clearvars info
    else
    fprintf([fn, '.h5 alrady exisits, skipping'])
    end
    
end

end