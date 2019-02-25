root=fullfile('/oak','stanford','groups','giocomo','attialex','NP_DATA');
Files= dir(fullfile(root,'*.mat'));

for iF=1:length(Files)
    try
    fprintf('Now working on %s \n',Files(iF).name )
    load(fullfile(root,Files(iF).name));
    good_cells = sp.cids(sp.cgs==2);
    idx=ismember(sp.clu,good_cells);

    spikes = double(sp.clu(idx));
    tempSP=[ones(size(spikes)) spikes spikes];
    mono=bz_MonoSynConvClick(double(tempSP),sp.st(idx),'plot',false);
    fprintf('Done for %s, now saving. \n',Files(iF).name)
    save(fullfile(root,Files(iF).name),'mono','-append','-v7.3')
    catch ME
    fprintf('Failed for %s \n',Files(iF).name)
    warning(ME.identifier)
    warning(ME.message)
    end
end

