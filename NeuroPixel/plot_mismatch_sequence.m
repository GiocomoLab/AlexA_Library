speed=true_speed;
speed_t=0.05;
% figure; plot(speed)
% 
all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
run_periods=smooth(speed,25)>speed_t;
run_window=-25:25;
valid_idx=true(size(mm_trigs));
for sort_idx=1:length(mm_trigs)
    runbi=sum(run_periods(mm_trigs(sort_idx)+run_window));
    if runbi<length(run_window)
        valid_idx(sort_idx)=false;
    end
end
% hold on
% plot(sit_periods)


mm_trigs=all_mm_trigs(valid_idx);
%%
[spike_mat,win,adata]=extract_triggered_spikes(sp,post(mm_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
%%
spike_mat=spike_mat(sp.cids+1,:,:);
%% turn spike_mat into firing rate
kernel=reshape(gausswin(401),1,1,[]);
zn=convn(spike_mat,kernel,'same');
%% separate into matrix for plotting and sorting, plus use only clusters with certain ID
cellIDX=sp.cgs>=1;
%cellIDX=true(1,size(zn,1));
mean_fr_sort=squeeze(mean(zn(cellIDX,2:2:end,:),2));
mean_fr_plot=squeeze(mean(zn(cellIDX,1:2:end,:),2));
%% plot and sort cells
[~,max_c]=max(mean_fr_sort,[],2);
mm_resp=mean(mean_fr_sort(:,4100:4599),2)-mean(mean_fr_sort(:,3500:3999),2);
[ss,sort_idx]=sort(mm_resp,'desc');
%[ss,ii]=sort(max_c);


[a,~]=max(mean_fr_plot,[],2);
ff_plot=bsxfun(@rdivide,mean_fr_plot,a);
ff_sort=bsxfun(@rdivide,mean_fr_sort,max(mean_fr_sort,[],2));

%figure('Name',filenames{iF})
figure
subplot(3,1,1)
imagesc(ff_sort(sort_idx,:))
title('data used for sorting')
subplot(3,1,2)
imagesc(ff_plot(sort_idx,:))
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


%%
cellID=sort_idx(1);
binned_array=squeeze(spike_mat(cellID,:,:));
bins=win(1):0.001:win(2);
    
    % set the new rasters
    [tr,b] = find(binned_array);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy*1+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
figure;subplot(2,1,1);plot(rasterX,rasterY)
hold on

subplot(2,1,2)
plot(-4:0.001:4,resp(cellID,:));
hold on
%%
figure
plot(squeeze(mean(adata,2)))
speed_diff=mean(adata(:,:,225:275),3)-mean(adata(:,:,150:200),3);
hold on
[a,b]=sort(speed_diff);
plot(squeeze(mean(adata(:,b(1:20),:),2)))
plot(squeeze(mean(adata(:,b(end-20:end),:),2)))

figure
resp = squeeze(mean(zn,2));
plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999)
resp = squeeze(mean(zn(:,b(1:20),:),2));
plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',[1 0 0])
resp = squeeze(mean(zn(:,b(end-20:end),:),2));
plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',[.5 1 0])

%%
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
depth=zeros(1,size(zn,1));
for iC=1:size(zn,1)
    depth(iC)=mean(spikeDepths(sp.clu==sp.cids(iC)));
end

%%
mean_fr_plot=squeeze(mean(zn,2));
[a,~]=max(mean_fr_plot,[],2);
ff_plot=bsxfun(@rdivide,mean_fr_plot,a);
[ss,sort_idx]=sort(depth,'desc');
figure;imagesc(ff_sort(sort_idx,:))
%%

figure
plot(speed)
hold on
temp=squeeze(mean(adata,2));
kernel=temp(150:250);
speed_matched=conv(speed,kernel,'same');
figure
plot(speed_matched);hold on;
%%
[pks,locs]=findpeaks(speed_matched,'MinPeakProminence',4);
figure
plot(speed_matched)
hold on
plot(locs,speed_matched(locs),'ro');
%plot(all_mm_trigs,speed_matched(all_mm_trigs),'ro')
trigs=locs+50;
[spike_mat2,win,adata2]=extract_triggered_spikes(sp,post(trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
%%
kernel=reshape(gausswin(401),1,1,[]);
zn2=convn(spike_mat2,kernel,'same');
figure
resp = squeeze(mean(zn2,2));
plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999)
