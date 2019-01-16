trials=[14:85];
spC=[];
dwell_time=[];
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    edges=[0:2:402];
    edges(1)=-5;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges);
    spC=cat(3,spC,spM);
    dwell_time=cat(1,dwell_time,dT);
end

spC=spC(sp.cgs>=1,:,:);
gauss_filter = gausswin(11)/sum(gausswin(11));

%%
avgMM=aggregateData.avgMM(aggregateData.AID==5,:);
mmresp=mean(avgMM(:,205:255),2)-mean(avgMM(:,141:190),2);
[a,b]=sort(mmresp,'desc');

if length(b)~= size(spC,1)
    error('Sure you have matching MM and spatial data??')
end
%%
for iC=1:24
    idx=b(iC);
fr=squeeze(spC(idx,:,:))'/0.02;
fr=fr./dwell_time;

figure('Position',[680    99   560   879]); subplot(5,1,1:2), 
imagesc(squeeze(fr))

%xlabel(sprintf('n spikes %d',sum(sum(spC(idx,:,:)))))

subplot(5,1,3)
fr=mean(spC(idx,:,:),3)/0.02;
dwell_t=mean(dwell_time,1);

fr_smoothed=conv(fr./dwell_t,gauss_filter,'same');
plot(0:200,fr_smoothed)
xlabel('averaged spatial fr')

subplot(5,1,4)
plot(-4:0.02:4,avgMM(idx,:))
xlabel('MM response')
subplot(5,1,5)
cid=idx-1;plot(sp.st(sp.clu==cid),sp.tempScalingAmps(sp.clu==cid),'.')
xlabel('spike amp during gain trials')
end
%%
