trials=[1:max(trial)];
%trials = trials(trial_gain == 1 & trial_contrast == 100);
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
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.cids)+1);
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

% 
% y=pdist(squeeze(spC(2,:,:))','spearman');
% Z1=linkage(y,'average');
% order=optimalleaforder(Z1,y);
% figure
% imagesc(ff(order,order));

%%
correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
diagAll=zeros(size(spatialMap,1),size(spatialMap,3)-1);
for iC=1:size(spatialMap,1)
    tmp=corr(squeeze(spatialMap(iC,:,:)));
    correlation_All(:,:,iC)=tmp;
    diagAll(iC,:)=diag(tmp,1);
end
diagAll(isnan(diagAll))=0;
figure
subplot(2,1,1)
imagesc(diagAll);
title('All 1st diags')
[coeff,score,latent,~,vex,mu] = pca(diagAll);
load1=score(:,1);
[a,b]=sort(load1);
subplot(2,1,2)
imagesc(diagAll(b,:))
title('Sorted by 1st pc')
print(gcf,fullfile('C:','tmp',[Files(iF).name(1:end-4) '_pca.png']),'-dpng')

figure
pp=nanmean(correlation_All,3);
imagesc(nanmean(correlation_All,3),[-.1 0.3])
title('mean over all corr')
transitions=find(diag(pp,1)<.10);
hold on
plot(transitions+1,transitions,'ro')
%% plot correlation map with trial structure image
trial_color = zeros(1,max(trial),3);
baseline = [0,0,0];
contrast = [.4 .4 .4];
%TODO extend for all different gains
gain_1 = [0 1 1];
gain_2 = [1 0 1];
for ii = 1:size(trial_color,2)
    if exist('trial_contrast','var')
        if trial_contrast(ii)<100
            trial_color(1,ii,:)=contrast;
        end
    end
    if exist('trial_gain','var')    
        if trial_gain(ii)<1
            trial_color(1,ii,:)=gain_1;
        end
        if trial_gain(ii)<.7
            trial_color(1,ii,:)=gain_2;
        end
    end
end
tc=permute(trial_color,[2 1 3]);
figure
subplot('Position',[0.1 0.1 .8 .8 ])
imagesc(nanmean(correlation_All,3),[-.1 0.5])
    xlabel(Files(iF).name,'Interpreter','None')

subplot('Position',[0.1 0.9 .8 0.05 ])
imagesc(trial_color)
set(gca,'XTick',[],'YTick',[])
subplot('Position',[0.05 0.1 .05 0.8 ])
imagesc(tc)
set(gca,'XTick',[],'YTick',[])
%print(gcf,fullfile('C:','tmp',[Files(iF).name(1:end-4) '_heatmap.png']),'-dpng')
%ca
