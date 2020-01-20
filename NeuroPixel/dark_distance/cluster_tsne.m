files = dir('F:\temp/data/np*dark*');

ACGAll=[];

for iF=1:numel(files)
    load(fullfile(files(iF).folder,files(iF).name))
    ACGAll=cat(1,ACGAll,ACG);
end

%%

ff=tsne(ACGAll);
%%
idx = ~isnan(ACGAll(:,1));
figure
subplot(1,2,1)
scatter3(ff(:,1),ff(:,2),log(PEAKValue(idx)./QUANT(idx)),3,PEAKS(idx))
zlim([-10 10])
subplot(1,2,2)
scatter(ff(:,1),ff(:,2),3,PXXPeak(idx))

%%
figure
while(true)
    subplot(1,2,1)
scatter(ff(:,1),ff(:,2),'.')
[xx,yy]=ginput()


sel_idx = ff(:,1)>xx(1) & ff(:,1)<xx(2) & ff(:,2)>yy(1) & ff(:,2)<yy(2);

hold on
scatter(ff(sel_idx,1),ff(sel_idx,2),'r.')

subplot(1,2,2)
plot(mean(ACGAll(sel_idx,:)))
pause
clf
end
%%

[cluidx]=kmeans(ACGAll(~isnan(ACGAll(:,1)),:),5);

%%
figure
scatter(ff(:,1),ff(:,2),3,cluidx)
%%
figure
tmp = ACGAll(~isnan(ACGAll(:,1)),:);
for ii=1:5
    hold on
    plot(mean(tmp(cluidx==ii,:)))
end