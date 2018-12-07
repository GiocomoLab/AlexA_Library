trials=[1:44];
spC=[];
dwell_time=[];
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',[0:2:402]);
    spC=cat(3,spC,spM);
    dwell_time=cat(1,dwell_time,dT);
end

% figure
% subplot(2,1,1)
% imagesc(spM)
% subplot(2,1,2)
% plot(dwell_time)

%%

smoothSigma = params.SmoothSigmaFR/params.SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

%% get spatial firing rate
iC=2;
fr=mean(spC(iC,:,:),3)/0.02;
dwell_t=mean(dwell_time,1);

fr_smoothed=conv(fr./dwell_t,gauss_filter,'same');
figure
hold on
plot(fr_smoothed)


%%
for iC=1:12
fr=squeeze(spC(iC,:,:))'/0.02;
fr=fr./dwell_time;

figure;imagesc(squeeze(fr))
xlabel(sprintf('n spikes %d',sum(sum(spC(iC,:,:)))))
end
%%
mean_fr_sort=mean(spC(good_cells,:,1:2:end),3);
avgDT_sort=mean(dwell_time(1:2:end,:));
mean_fr_sort=mean_fr_sort./avgDT_sort/0.02;
mean_fr_sort=conv2(mean_fr_sort,gauss_filter','same');



mean_fr_plot = mean(spC(good_cells,:,2:2:end),3);
avgDT_plot = mean(dwell_time(2:2:end,:));
mean_fr_plot = mean_fr_plot./avgDT_plot/0.02;
mean_fr_plot = conv2(mean_fr_plot,gauss_filter','same');




[~,max_c]=max(mean_fr_sort,[],2);

[ss,ii]=sort(max_c);

[a,~]=max(mean_fr_plot,[],2);
ff_plot=bsxfun(@rdivide,mean_fr_plot,a);
ff_sort=bsxfun(@rdivide,mean_fr_sort,max(mean_fr_sort,[],2));

figure('Name','bla')
subplot(3,1,1)
imagesc(ff_sort(ii,:))
title('data used for sorting')
subplot(3,1,2)
imagesc(ff_plot(ii,:))
title('held out data')
