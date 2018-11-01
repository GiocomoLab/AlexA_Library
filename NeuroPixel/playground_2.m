
speed = diff(posx);
%figure;plot(posx)
jumps=find(speed<-40);
for ii=1:length(jumps)
    speed(jumps(ii))=.5*speed(jumps(ii)-1)+.5*speed(jumps(ii)+1);
end
speed_t=0.1;
figure; plot(speed)

sit_periods=speed<speed_t;
hold on
plot(sit_periods)

%%
transitions=find(diff([0;sit_periods;0]));
sit_times=transitions(2:2:end)-transitions(1:2:end);

figure
histogram(sit_times,80)

set(gca,'XTick',[0:25:1600],'XTickLabel',[0:25:1600]/50)

sampling_rate=50;
sitwin=[zeros(1,4*sampling_rate) ones(1,2.5*sampling_rate)];
stop_idx=strfind(sit_periods',sitwin)+4*sampling_rate;
stop_time=post(stop_idx);
%%
[spike_mat,win,adata]=extract_triggered_spikes(sp,stop_time,'win',[-4 4],'aux',[post'; [0 speed']],'aux_win',[-200 200]);
%%

%%
ff=reshape(gausswin(81),1,1,[]);
zn=convn(spike_mat,ff,'same');
%%
mean_fr_sort=squeeze(mean(zn(:,2:2:end,:),2));
mean_fr_plot=squeeze(mean(zn(:,1:2:end,:),2));
%%
[~,max_c]=max(mean_fr_sort,[],2);

[ss,ii]=sort(max_c);

[a,~]=max(mean_fr_plot,[],2);

ff=bsxfun(@rdivide,mean_fr_plot,a);

figure
imagesc(ff(ii,:))
