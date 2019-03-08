root=fullfile('/oak','stanford','groups','giocomo','attialex','NP_DATA');
Files= dir(fullfile(root,'*.mat'));

for iF=1:length(Files)
    load(fullfile(root,Files(iF).name),'connected');
    try
    fprintf('Now working on %s \n',Files(iF).name )
    load(fullfile(root,Files(iF).name));
    glm_playground
    fprintf('Done for %s, now saving. \n',Files(iF).name)
    save(fullfile(root,Files(iF).name),'glmModel','-append')
    catch ME
    fprintf('Failed for %s \n',Files(iF).name)
    warning(ME.message)
    end
end
