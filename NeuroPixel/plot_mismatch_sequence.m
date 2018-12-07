speed=true_speed;
speed_t=0.05;
% figure; plot(speed)
% 
mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
run_periods=smooth(speed,25)>speed_t;
run_window=-25:25;
valid_idx=true(size(mm_trigs));
for ii=1:length(mm_trigs)
    runbi=sum(run_periods(mm_trigs(ii)+run_window));
    if runbi<length(run_window)
        valid_idx(ii)=false;
    end
end
% hold on
% plot(sit_periods)


mm_trigs=mm_trigs(valid_idx);
%%
[spike_mat,win,adata]=extract_triggered_spikes(sp,post(mm_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
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
%%
resp = squeeze(mean(zn,2));
%%
params=struct();
params.winIDX=-4000:4000;
params.masterTime=params.winIDX/1000;
params.xLim=[-1 3];
figure

plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999)
figure


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
