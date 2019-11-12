chunksize=50;
stride = 20;
files = dir('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*_gain*.mat');
trials = [11:34];
PEAKS=zeros(24,8,numel(files));
SHIFTS = PEAKS;
%parpool(8)
parfor iF = 5:numel(files)
    try
    data = load(fullfile(files(iF).folder,files(iF).name));
    [peak,shift]=calculatePeakShiftSession(data,trials,chunksize,stride,'MEC',0.2);
    PEAKS(:,:,iF)=peak;
    SHIFTS(:,:,iF)=shift;
    startVec = 1:stride:(200-chunksize);
    
    x=startVec+25;
    x = x-1;
    x = x*2;
    fig = figure('visible','off');
    hold on
for iT = 1:24
    tmp = x+400*(iT-1);
    if ismember(iT,[11:14])
            plot(tmp,peak(iT,:),'r.')

    else
        
    plot(tmp,peak(iT,:),'b.')
    end
end
    ylim([0.2, 0.8])
    drawnow
    session_name = files(iF).name(1:end-4);
    saveas(fig,fullfile('/oak/stanford/groups/giocomo/attialex/Images/xcorrv1',[session_name +'.png']))
    close(fig)
    catch ME
       disp(ME.message)
       session_name = files(iF).name(1:end-4);
       fprintf('not working for %s \n',session_name)
    end
    
end