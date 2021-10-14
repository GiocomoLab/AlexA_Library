savepath = 'F:\Alex\mm_out\AA_2*';


savefiles=dir(savepath);
plot_path = fullfile(savefiles(1).folder,'rasters');
if ~isfolder(plot_path)
    mkdir(plot_path);
end

for iF=1:numel(savefiles)
    mf = matfile(fullfile(savefiles(iF).folder,savefiles(iF).name));
    if nnz(strcmp(mf.region,'ENTm'))==0
        continue
    end
    plot_root = fullfile(plot_path,savefiles(iF).name(1:end-4));
    if ~isfolder(plot_root)
        mkdir(fullfile(plot_root));
    end
    spike_times = mf.spike_times;
    region = mf.region;
    vec = cat(1,spike_times{:});
    gc = unique(vec(:,2));
    fig=figure('Position',[680   851   560   127],'Visible','Off');
    for iC=1:numel(gc)
        c_idx = vec(:,2)==gc(iC);
        scatter(vec(c_idx,1),vec(c_idx,3),15,'k','.')
        xline(0)
        savename = fullfile(plot_root,sprintf('%s_%d.png',region{iC},gc(iC)));
        saveas(gcf,savename)
        cla
    end
    close(fig)
end