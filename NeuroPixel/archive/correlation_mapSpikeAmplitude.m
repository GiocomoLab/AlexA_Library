trials=[1:max(trial)];

spatialMap=[];
dwell_time=[];
edges=[0:10:410];
edges(1)=-.01;
posx(posx<0)=0;
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.clu)+1);
    spatialMap=cat(3,spatialMap,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
spatialMap=spatialMap(sp.cids+1,:,:);
spatialMap=spatialMap(sp.cgs==2,:,:);
spatialMap=spatialMap(:,1:end-1,:);
dwell_time=dwell_time(:,1:end-1);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end

correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
diagAll=zeros(size(spatialMap,1),size(spatialMap,3)-1);
for iC=1:size(spatialMap,1)
    tmp=corr(squeeze(spatialMap(iC,:,:)));
    correlation_All(:,:,iC)=tmp;
end

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
cells_to_plot = sp.cids(sp.cgs==2);
amplitudes = zeros(length(cells_to_plot),max(trial));
for iC=1:length(cells_to_plot)
spike_id=sp.clu==cells_to_plot(iC);
spike_t = sp.st(spike_id);
spike_idx = discretize(spike_t,post);

for iT=1:max(trial)
    tmp_idx = spike_idx & trial(spike_idx) == iT;
    amplitudes(iC,iT)=mean(spikeAmps(tmp_idx));
end
end


tmp=nanmean(correlation_All,3);
imagesc(tmp)
set(gca,'CLim',prctile(tmp(:),[20 80]))
yyaxis right
tmpA=nanmean(amplitudes);
plot(tmpA,'k','LineWidth',2)
ylim([min(tmpA) max(tmpA)+1*range(tmpA)]);
