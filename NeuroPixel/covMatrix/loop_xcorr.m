chunksize=50;
stride = 20;
files = dir('F:\NP_DATA\npF4*_gain*.mat');
trials = [11:34];
PEAKS=zeros(24,8,numel(files));
SHIFTS = PEAKS;
for iF = 1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    [peak,shift]=calculatePeakShiftSession(data,trials,chunksize,stride,'MEC',0.2);
    PEAKS(:,:,iF)=peak;
    SHIFTS(:,:,iF)=shift;
    startVec = 1:stride:(200-chunksize);
    
    x=startVec+25;
    x = x-1;
    x = x*2;
    fig = figure();
    hold on
for iT = 1:24
    tmp = x+400*(iT-1);
    plot(tmp,peak(iT,:),'b.')
end
    ylim([0, 1]])
    drawnow
    session_name = files(iF).name(1:end-4);
    %saveas(fig,fullfile('/oak/stanford/groups/giocomo/attialex/Images/xcorrv1',[session_name +'.png']))
    %close(fig)
end