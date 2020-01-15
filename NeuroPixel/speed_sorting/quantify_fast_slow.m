files = dir('Z:\giocomo\attialex\speed_sort4\*.mat');
delays = [];
region = {};
slow_trials = [];
fast_trials = [];
for iF=1:numel(files)
    load(fullfile(files(iF).folder,files(iF).name));
    delays = cat(2,delays,data.delay);
    region = cat(2,region,data.region);
    slow_trials = cat(2,slow_trials,data.slow_trials);
    fast_trials = cat(2,fast_trials,data.fast_trials);

end

%%
[aa,max_i]=max((fast_trials+slow_trials)/2);

reg_idx = startsWith(region,'VISp') & delays(2,:)>0.01;

figure
n_B=4;
step = 50;
for ii = 1:n_B
subplot(1,n_B,ii)
this_idx = reg_idx & (ii-1)*step+1<max_i & max_i<ii*step;
tmp = mean(slow_trials(:,this_idx),2);
tmp = tmp/max(tmp);
plot(tmp)
hold on
tmp = mean(fast_trials(:,this_idx),2);
tmp = tmp/max(tmp);
plot(tmp)
title(mean(delays(1,this_idx)))
end
%%
delays(1,:)=round(delays(1,:));
%%
idx = startsWith(region,'VISp') & delays(2,:)>0.01;
figure
h1 = histogram(delays(1,idx));

hold on
idx = startsWith(region,'MEC') & delays(2,:)>0.01;
h2 = histogram(delays(1,idx));

h1.Normalization = 'probability';
h1.BinLimits=[-10.5 10.5];
h1.NumBins=21;
h2.Normalization = 'probability';
h2.BinLimits=[-10.5 10.5];
h2.NumBins=21;

legend({'V1','MEC'})

%%

avg = (slow_trials+fast_trials/2);

avgNorm = avg-mean(avg);
idx = startsWith(region,'CA') & delays(2,:)>0.01;

visTC=avgNorm(:,idx)';
nClu=6;
clu = kmeans(visTC,nClu,'Replicates',5);
figure
hold on
slow_vc=slow_trials(:,idx);
fast_vc=fast_trials(:,idx);
for ii=1:nClu
    plot(nanmean(visTC(clu==ii,:)))
end
figure

delays_this = delays(1,idx);
for ii=1:nClu
    subplot(1,nClu,ii)
    iidx = clu==ii;
    plot(nanmean(slow_vc(:,iidx),2))
    hold on
    plot(nanmean(fast_vc(:,iidx),2))
    title(mean(delays_this(iidx)));
end
    