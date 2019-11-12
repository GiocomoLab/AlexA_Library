chunksize=50;
stride = 20;
files = dir('F:\NP_DATA\npF4*_gain*.mat');
for iF = 1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    [PEAKS,SHIFTS]=calculatePeakShiftSession(data,trials,chunksize,stride);
     x=startVec+25;
    x = x-1;
    x = x*2;
    figure
    hold on
for iT = 1:24
    tmp = x+400*(iT-1);
    plot(tmp,PEAKS(iT,:),'b.')
end
    saveas(fig,fullfile('/oak/stanford/groups/giocomo/attialex/Images/xcorrv1',[session_name +'.png']))
end