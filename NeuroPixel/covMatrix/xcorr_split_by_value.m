gains = {'1.0','0.8','0.5'};
for ii=2:3%1:3
    load(['Z:\giocomo\attialex\Images\xcorrv7\allData_MEC_' gains{ii} '_100.mat'])
nans = isnan(output.Y(:,1));
validY=output.Y(~nans,:);

meanY=mean(validY(:,61:70),2);
[a,sid]=sort(meanY);

zz=zeros(16*200,size(output.SPEED,3));
for iT=1:16
    tmp = squeeze(output.FR(iT,:,:));
    if ismember(iT,[7:10])
        %tmp=tmp/str2double(gains{ii});
    end
    idx = ((iT-1)*200)+1:iT*200;
    zz(idx,:)=tmp;
end
zz=bsxfun(@rdivide,zz,nanmean(zz(1:6*200,:)));
xvec = (2:2:(16*400))-6*400;

figure

nSlice = 2;
nperSlice = floor(size(validY,1)/nSlice);
for iS=1:nSlice
    subplot(2,1,1)
    hold on
    idx = ((iS-1)*nperSlice+1):iS*nperSlice;
    plot(output.X-6*400,mean(validY(sid(idx),:)),'-')
    subplot(2,1,2)
    hold on
    plot(xvec,squeeze(nanmean(zz(:,idx),2)));
end
subplot(2,1,1)
ylim([0 0.6])
xlabel('Distance from onset [cm]')
ylabel('Peak xcorr')
title(['Gain = ' gains{ii}])
subplot(2,1,2)
xlim([-200 1200])
end