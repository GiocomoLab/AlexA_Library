 load('/oak/stanford/groups/giocomo/attialex/MM_aggregate12-Sep-2019.mat')
 %%
 addpath(genpath('/home/users/attialex/AlexA_Library'))
 addpath(genpath('/home/users/attialex/boundedline-pkg'))
 %% plot visual cortex
IDX = strcmp(aggregateData.PARENT,'VISp')' & aggregateData.CGS>1;
params=struct();
params.winIDX=-200:200;
params.masterTime=params.winIDX/50;
params.xLim=[-1 3];
figure

plotAVGSEM(aggregateData.avgMM(IDX,:)',gca,'parameters',params,'ms',true,'baseline',165:190)
xlabel('time [s]')
ylabel('Firing rate change [Hz]')
grid on
% 
% [a,b]=sort(aggregateData.DEPTH);
% figure
% plotAVGSEM(aggregateData.avgMM(b(1:200),:)',gca,'parameters',params,'ms',true,'baseline',165:190)
% plotAVGSEM(aggregateData.avgMM(b(end-200:end),:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])
% xlabel('time [s]')
% ylabel('Firing rate change [Hz]')
% grid on
% legend({'ventral','Dorsal'})

%%
[a,b]=sort(aggregateData.DEPTH(IDX));
%ff_plot=bsxfun(@rdivide,aggregateData.avgMM(aggregateData.CGS==2,150:350),a);
ff_plot=aggregateData.avgMM-mean(aggregateData.avgMM(:,170:195),2);
ff_plot=ff_plot(IDX,:);
figure
imagesc(ff_plot(b,150:350),[-10 10])
set(gca,'XTick',[0:50:200])
set(gca,'XTickLabel',[-1:1:4])
colormap jet
colorbar
xlabel('Time, [s]')
ylabel('Unit #')

%%
%avgMM = aggregateData.avgMM(IDX,:);

%%
IDXSite = aggregateData.AID==3;
IDX1= true(nnz(aggregateData.AID==3),1);

IDX = strcmp(aggregateData.PARENT(IDXSite),'VISp')' & aggregateData.CGS(IDXSite)>1;
avgMM = aggregateData.avgMM(IDXSite,:);
avgMM=avgMM(IDX,:);
mResp = mean(avgMM(:,205:250),2)-mean(avgMM(:,175:190),2);
[a,b]=sort(mResp);
ff_plot=avgMM-mean(avgMM(:,170:195),2);

figure
imagesc(ff_plot(b,150:350),[-10 10])
set(gca,'XTick',[0:50:200])
set(gca,'XTickLabel',[-1:1:4])
colormap jet
colorbar
xlabel('Time, [s]')
ylabel('Unit #')


MM_snps = aggregateData.MM_snps{3}(IDX,:,:);
for ii = 1:5
    figure
    idx=b(ii);
    imagesc(squeeze(MM_snps(idx,:,:)));
end