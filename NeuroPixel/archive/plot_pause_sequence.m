
speed = diff(posx);
%figure;plot(posx)
jumps=find(speed<-40);
for ii=1:length(jumps)
    speed(jumps(ii))=.5*speed(jumps(ii)-1)+.5*speed(jumps(ii)+1);
end
speed_t=0.2;
% figure; plot(speed)
% 
sit_periods=speed<speed_t;
% hold on
% plot(sit_periods)

%%
transitions=find(diff([0;sit_periods;0]));
sit_times=transitions(2:2:end)-transitions(1:2:end);

%figure
%histogram(sit_times,80)

%set(gca,'XTick',[0:25:1600],'XTickLabel',[0:25:1600]/50)

sampling_rate=50;
sitwin=[zeros(1,4*sampling_rate) ones(1,2.5*sampling_rate)];
stop_idx=strfind(sit_periods',sitwin)+4*sampling_rate;
stop_time=post(stop_idx);
%%
[spike_mat,win,adata]=extract_triggered_spikes(sp,stop_time,'win',[-4 4],'aux',[post'; [0 speed']],'aux_win',[-200 200]);
%%

%% turn spike_mat into firing rate
kernel=reshape(gausswin(401),1,1,[]);
zn=convn(spike_mat,kernel,'same');
%% separate into matrix for plotting and sorting, plus use only clusters with certain ID
cellIDX=sp.cids(sp.cgs==2)+1;
mean_fr_sort=squeeze(mean(zn(cellIDX,2:2:end,:),2));
mean_fr_plot=squeeze(mean(zn(cellIDX,1:2:end,:),2));
%% plot and sort cells
[~,max_c]=max(mean_fr_sort,[],2);

[ss,ii]=sort(max_c);

[a,~]=max(mean_fr_plot,[],2);
ff_plot=bsxfun(@rdivide,mean_fr_plot,a);
ff_sort=bsxfun(@rdivide,mean_fr_sort,max(mean_fr_sort,[],2));

%figure('Name',filenames{iF})
figure
subplot(3,1,1)
imagesc(ff_sort(ii,:))
title('data used for sorting')
subplot(3,1,2)
imagesc(ff_plot(ii,:))
title('held out data')
subplot(3,1,3)
plot(-4:0.02:4,squeeze(adata))
%% extract firing rate for each neuron, and correlate with velocity to get cells with running correlation < 0 
[sp_Mat,zn_all,zn_all_d]=get_spike_mat(sp);

aa=corr(zn_all_d',speed(1:end));

cellIDX=sp.cids(sp.cgs==2 & aa(sp.cids+1)'<0)+1;
mean_fr_sort=squeeze(mean(zn(cellIDX,2:2:end,:),2));
mean_fr_plot=squeeze(mean(zn(cellIDX,1:2:end,:),2));
%% plot and sort cells for cells with running corr < 0 
[~,max_c]=max(mean_fr_sort,[],2);

[ss,ii]=sort(max_c);

[a,~]=max(mean_fr_plot,[],2);

%%
ff_plot=bsxfun(@rdivide,mean_fr_plot,a);
ff_sort=bsxfun(@rdivide,mean_fr_sort,max(mean_fr_sort,[],2));

figure('Name',filenames{iF})
subplot(3,1,1)
imagesc(ff_sort(ii,:))
title('data used for sorting')
subplot(3,1,2)
imagesc(ff_plot(ii,:))
title('held out data')
subplot(3,1,3)
plot(-4:0.02:4,squeeze(adata))

%% get spike depth

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
depth=zeros(1,size(sp_Mat,1));
for iC=1:size(sp_Mat,1)
    depth(iC)=mean(spikeDepths(sp.clu==iC-1));
end
%%
figure('Name',filenames{iF})
plot(aa,depth,'.')
hold on
plot(aa(sp.cgs==2),depth(sp.cgs==2),'r.')
xlabel('correlation with speed')
ylabel('depth on probe')
legend({'MUA','Good'})
