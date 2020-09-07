speed = diff(posx);
%figure;plot(posx)
jumps=find(speed<-40);
for ii=1:length(jumps)
    speed(jumps(ii))=.5*speed(jumps(ii)-1)+.5*speed(jumps(ii)+1);
end
%%
ff=abs(posx-250);

ff(ff>.5)=2;
ff(ff<.5)=0;

TF=islocalmin(ff,'FlatSelection', 'center');

figure;plot(post,posx,post(TF),posx(TF),'r*')
stop_time=post(TF);
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

figure('Name',filenames{iF})
subplot(3,1,1)
imagesc(ff_sort(ii,:))
title('data used for sorting')
subplot(3,1,2)
imagesc(ff_plot(ii,:))
title('held out data')
subplot(3,1,3)
plot(-4:0.02:4,squeeze(adata))