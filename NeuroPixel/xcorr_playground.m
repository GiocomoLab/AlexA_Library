trials=[1:max(trials)];
spC=[];
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
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges);
    spC=cat(3,spC,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
spC=spC(sp.cids+1,:,:);
spC=spC(sp.cgs==2,:,:);
spC=spC(:,1:end-1,:);
dwell_time=dwell_time(:,1:end-1);
%%
gauss_filter = gausswin(5)/sum(gausswin(5));

dt=dwell_time';

dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spC,1)
    spC(ii,:,:)=spC(ii,:,:)./dt;
end
spC(isnan(spC))=0;
dia=diag(ff,1);

figure
imagesc(squeeze(spC(1,:,:))')
ff=corr(squeeze(spC(1,:,:)));
transitions=find(dia<0);
figure
imagesc(ff)
hold on
plot(transitions+1,transitions,'ro')
figure
plot(dia);

y=pdist(squeeze(spCd(1,:,:))','spearman');
Z1=linkage(y,'average');
order=optimalleaforder(Z1,y);
figure
imagesc(ff(order,order));

%%
ffAll=zeros(size(ff,1),size(ff,2),size(spC,1));
diagAll=zeros(size(spC,1),size(ff,1)-1);
for iC=1:size(spC,1)
    tmp=corr(squeeze(spC(iC,:,:)));
    ffAll(:,:,iC)=tmp;
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
figure
imagesc(nanmean(ffAll,3),[-.1 0.3])
title('mean over all corr')
%%
for i
