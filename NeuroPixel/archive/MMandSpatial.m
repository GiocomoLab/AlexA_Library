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
%cellIDX=find(sp.cgs>=1);
spC=spC(sp.cids+1,:,:);
%%
dt=dwell_time';

dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spC,1);
    spC(ii,:,:)=spC(ii,:,:)./dt;
end
    
gauss_filter = gausswin(11)/sum(gausswin(11));

%%
avgMM=aggregateData.avgMM(aggregateData.AID==2,:);
mmresp=mean(avgMM(:,205:255),2)-mean(avgMM(:,141:190),2);
[a,b]=sort(mmresp,'desc');

if length(b)~= size(spC,1)
    error('Sure you have matching MM and spatial data??')
end
%%
for iC=1:20
    idx=stabi(iC);
% fr=squeeze(spC(idx,:,:))'/0.02;
% fr=fr./dwell_time;
fr=squeeze(spC(idx,:,:));
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
%plot(-4:0.02:4,avgMM(idx,:))
params=struct();
params.winIDX=-200:200;
params.masterTime=params.winIDX/50;
params.xLim=[-1 3];

%plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',175:195)
plot(params.masterTime,avgMM(idx,:));
xlim([-1 3])
xlabel('MM response')
subplot(5,1,5)
% cid=sp.cids(idx);plot(sp.st(sp.clu==cid),sp.tempScalingAmps(sp.clu==cid),'.')
% xlabel('spike amp during gain trials')

binned_array=squeeze(spike_mat(idx,:,:));
bins=win(1):0.001:win(2);
    
    % set the new rasters
    [tr,b] = find(binned_array);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy*1+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
plot(rasterX,rasterY)
hold on

end
%%
