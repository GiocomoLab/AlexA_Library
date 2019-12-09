gains = {'1.0','0.8','0.5'};
for ii=1:2%1:3
    load(['Z:\giocomo\attialex\Images\xcorrv7\allData_MEC_' gains{ii} '_100.mat'])
nans = isnan(output.Y(:,1));
validY=output.Y(~nans,:);
validS = output.S(~nans,:);
T = 10;
meanY=mean(validY(:,(T-1)*10+1:T*10),2);
[a,sid]=sort(meanY);

run_speed=zeros(16*200,nnz(~nans));
firing_rate = run_speed;
for iT=1:16
    tmp = squeeze(output.SPEED(iT,:,~nans));
    tmpfr = squeeze(output.FR(iT,:,~nans));
    if ismember(iT,[7:10])
        tmp=tmp/str2double(gains{ii});
    end
    idx = ((iT-1)*200)+1:iT*200;
    run_speed(idx,:)=tmp;
    firing_rate(idx,:)=tmpfr;
end
run_speed=bsxfun(@rdivide,run_speed,nanmean(run_speed(1:6*200,:)));
firing_rate=bsxfun(@rdivide,firing_rate,nanmean(firing_rate(1:6*200,:)));
xvec = (2:2:(16*400))-6*400;

figure

nSlice = 3;
nperSlice = floor(size(validY,1)/nSlice);
for iS=1:nSlice
    subplot(3,2,1)
    hold on
    idx = ((iS-1)*nperSlice+1):iS*nperSlice;
    plot(output.X-6*400,mean(validY(sid(idx),:)),'-')
    
    subplot(3,2,2)
    hold on
    plot(output.X-6*400,mean(validS(sid(idx),:)),'-')
    
    subplot(3,2,[3 4])
    hold on
    plot(xvec,squeeze(nanmean(run_speed(:,idx),2)));
    subplot(3,2,[5 6])
    hold on
    plot(xvec,squeeze(nanmean(firing_rate(:,idx),2)));
end
subplot(3,2,1)
ylim([0 0.6])
xlabel('Distance from onset [cm]')
ylabel('Peak xcorr')
title(['Gain = ' gains{ii}])
subplot(3,2,[3 4])
title('run speed')
xlim([-200 1200])
subplot(3,2,[5 6])
title('firing rate')
xlim([-200 1200])
end
%%

gains = {'1.0','0.8'};
for ii=1:numel(gains)
    load(['Z:\giocomo\attialex\Images\xcorrv7\allData_MEC_' gains{ii} '_100.mat'])
nans = isnan(output.Y(:,1));
validY=output.Y(~nans,:);
validS = output.S(~nans,:);
T2 = 8;
T1=7;
meanY=mean(validY(:,(T1-1)*10+1:T1*10-5),2)-mean(validY(:,(T1-1)*10+6:T1*10),2);
[a,sid]=sort(meanY);

run_speed=zeros(16*200,nnz(~nans));
firing_rate = run_speed;
for iT=1:16
    tmp = squeeze(output.SPEED(iT,:,~nans));
    tmpfr = squeeze(output.FR(iT,:,~nans));
    if ismember(iT,[7:10])
        tmp=tmp/str2double(gains{ii});
    end
    idx = ((iT-1)*200)+1:iT*200;
    run_speed(idx,:)=tmp;
    firing_rate(idx,:)=tmpfr;
end
run_speed=bsxfun(@rdivide,run_speed,nanmean(run_speed(1:6*200,:)));
firing_rate=bsxfun(@rdivide,firing_rate,nanmean(firing_rate(1:6*200,:)));
xvec = (2:2:(16*400))-6*400;

figure

nSlice = 3;
nperSlice = floor(size(validY,1)/nSlice);
for iS=1:nSlice
    subplot(3,2,1)
    hold on
    idx = ((iS-1)*nperSlice+1):iS*nperSlice;
    plot(output.X-6*400,mean(validY(sid(idx),:)),'-')
    
    subplot(3,2,2)
    hold on
    plot(output.X-6*400,mean(validS(sid(idx),:)),'-')
    
    subplot(3,2,[3 4])
    hold on
    plot(xvec,squeeze(nanmean(run_speed(:,idx),2)));
    subplot(3,2,[5 6])
    hold on
    plot(xvec,squeeze(nanmean(firing_rate(:,idx),2)));
end
subplot(3,2,1)
ylim([0 0.6])
xlabel('Distance from onset [cm]')
ylabel('Peak xcorr')
title(['Gain = ' gains{ii}])
subplot(3,2,2)
xlabel('Distance from onset [cm]')
ylabel('Shift [cm]')
title(['Shift'])

subplot(3,2,[3 4])
title('run speed')
xlim([-200 1200])
subplot(3,2,[5 6])
title('firing rate')
xlim([-200 1200])
end

%%
figure
for ii=1:12
    subplot(3,4,ii)
    %plot(meanY,validS(:,68+ii),'.')
    mdl = fitlm(meanY,validS(:,70+ii));
    plot(mdl)
    pp=anova(mdl);
    legend('off')
    xlabel(sprintf('%.3f',pp{1,5}))
    ylim([-10 10])
end