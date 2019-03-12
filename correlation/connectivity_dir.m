root=fullfile('/oak','stanford','groups','giocomo','attialex','NP_DATA');
Files= dir(fullfile(root,'*.mat'));

for iF=1:length(Files)
    clear connected
    load(fullfile(root,Files(iF).name),'connected');
    %if ~exist('connected','var')
    try
        fprintf('Now working on %s \n',Files(iF).name )
        load(fullfile(root,Files(iF).name));
        good_cells = sp.cids(sp.cgs==2);
        idx=ismember(sp.clu,good_cells);
        
        spikes = double(sp.clu(idx))+1; %add 1 bc it fails for cluID ==0
        tempSP=[ones(size(spikes)) spikes spikes];
        mono=bz_MonoSynConvClick(double(tempSP),sp.st(idx),'plot',false);
        connected = mono.sig_con;
        connected = connected -1;%subtract 1 to make connected correspond to cluid
        fprintf('Done for %s, now saving. \n',Files(iF).name)
        save(fullfile(root,Files(iF).name),'connected','-append')
    catch ME
        fprintf('Failed for %s \n',Files(iF).name)
        warning(ME.message)
    end
    %end
end

