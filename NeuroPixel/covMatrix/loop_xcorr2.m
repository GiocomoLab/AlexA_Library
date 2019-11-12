chunksize=50;
stride = 20;
trials = [11:34];

region = 'MEC';
contrast = 100;
gain_to_look_at = 0.8;

[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'F:/NP_DATA');
%%
n_chunks = 0;
for ii=1:numel(filenames)
    n_chunks = n_chunks+numel(triggers{ii});
end
PEAKS=zeros(24,8,n_chunks);
SHIFTS = PEAKS;
tt=(-10:13);
cntr = 0;
for iF = 1:numel(filenames)
    data = load(filenames{iF});
    for iTrigger = 1:numel(triggers{iF})
        cntr = cntr+1;
        current_trig = triggers{iF}(iTrigger);
        trials = current_trig+tt;
        
        [peak,shift]=calculatePeakShiftSession(data,trials,chunksize,stride,region,0.2);
        PEAKS(:,:,cntr)=peak;
        SHIFTS(:,:,cntr)=shift;
    startVec = 1:stride:(200-chunksize);
    
    x=startVec+25;
    x = x-1;
    x = x*2;
    fig = figure('visibile','off');
    hold on
for iT = 1:24
    tmp = x+400*(iT-1);
    plot(tmp,peak(iT,:),'b.')
end
    ylim([0.2, 0.8])
    drawnow
    [~,session_name,~] = fileparts(filenames{iF});
    savepath = '/oak/stanford/groups/giocomo/attialex/Images/xcorrv1';
    saveas(fig,fullfile(savepath,sprintf('%s_%s_%.1f_%d.png',session_name,region,gain_to_look_at,contrast)))
    close(fig)
    end
end
%%
X=[];

startVec = 1:stride:(200-chunksize);
    
    x=startVec+25;
    x = x-1;
    x = x*2;
for iT = 1:24
    tmp = x+400*(iT-1);
    X=cat(2,X,tmp);
    %plot(tmp,peak(iT,:),'b.')
end
Y=zeros(size(PEAKS,3),numel(X));
chunks_per_trial = size(PEAKS,2);
for iS = 1:size(PEAKS,3)
    for iT = 1:24
        idx = ((iT-1)*chunks_per_trial+1):iT*chunks_per_trial;
        tmp = squeeze(PEAKS(iT,:,iS));
        Y(iS,idx)=tmp;
        
    end
end