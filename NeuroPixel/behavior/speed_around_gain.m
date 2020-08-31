

ops.BinWidth =2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.TimeBin = 0.02;

fi = gausswin(11);
fi=fi'/sum(fi);
ops.filter = fi;
spfi = gausswin(5)'/sum(gausswin(5));
ops.speed_filter = spfi;
ops.n_trials = 10;
ops.plotfig = false;
ops.maxLag = 20; % in cm
%OAK='/oak/stanford/groups/giocomo/';
OAK='/Volumes/Samsung_T5/';
%% savedir =
savedir = fullfile(OAK,'attialex','speed_filtered_greedy');
%savedir = fullfile('F:/temp/','speed_filtered');
imdir = fullfile(savedir,'images');



%% find files


gain = 0.5;
contrast = 100;
region = 'VISp';
regions = {'VISp','MEC','RS'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
%%
nTot=0;
for iF=1:numel(filenames)
    for iT =1:numel(triggers{iF})
        nTot = nTot+1;
    end
end
%%
trials2show = [-6:9];
nT=numel(trials2show);
speedMat = nan(nT,ops.nBins,nTot);
speedMatRaw = nan(nT,ops.nBins,nTot);
cntr=0;
for iF=1:numel(filenames)
    data = load(filenames{iF},'posx','post','trial_gain','trial');
    [speed,speed_raw]=calcSpeed(data.posx,ops);
    %nT=numel(data.trial_gain);
    
    
    
    
    
    discrete_pos = discretize(data.posx,ops.edges);
    
    for iTrigger=1:numel(triggers{iF})
        cntr = cntr+1;
        trials = triggers{iF}(iTrigger)+trials2show;
        trial_cntr = 0;
        for iT=trials
            trial_cntr = trial_cntr+1;
            for iB=1:ops.nBins
                idx = data.trial==iT & discrete_pos == iB;
                if nnz(idx)>0
                    val = mean(speed(idx));
                    val_raw = mean(speed_raw(idx));
                    speedMat(trial_cntr,iB,cntr)=val/data.trial_gain(iT);
                    speedMatRaw(trial_cntr,iB,cntr)=val_raw/data.trial_gain(iT);
                end
            end
        end
    end
end
%%

% speedMat_c=speedMat;
% speedMat_c(7:10,:,:)=speedMat_c(7:10,:,:)/.8;
tmp=squeeze(nanmean(speedMat,1));
normfact =squeeze(nanmean(tmp,2));
normfact = max(tmp);
speedMatNorm = speedMat;
speedMatNormRaw = speedMatRaw;
for iRep = 1:size(speedMat,3)
    speedMatNorm(:,:,iRep)=speedMatNorm(:,:,iRep)./normfact(iR);
end
%%
figure
imagesc(nanmean(speedMatRaw,3))


%%
figure
hold on
 cols = cbrewer('qual','Dark2',16);
for iT = 1:size(speedMatNorm,1)
    tmp = nanmean(speedMat(iT,:,:),3);
    plot(tmp,'-','Color',cols(iT,:))
end
%%
t1=nanmean(speedMatNorm(6,:,:),3);
t2=nanmean(speedMatNorm(7,:,:),3);
figure
plot(t1,'k')
hold on
plot(t2,'Color',get_color(0.8,100))
  %%
  
  t1=nanmean(speedMatRaw(6,:,:),3);
t2=nanmean(speedMatRaw(7,:,:),3);
figure
plot(t1,'k')
hold on
plot(t2,'Color',get_color(0.5,100))
%%
  t1=nanmean(nanmean(speedMatRaw(3:6,:,:),3),1);
t2=nanmean(nanmean(speedMatRaw(7:10,:,:),3),1);
t3 = nanmean(nanmean(speedMatRaw(11:15,:,:),3),1);
figure
plot(t1,'k')
hold on
plot(t2,'Color',get_color(0.8,100))
plot(t3,'k--')
%%
tmpSM = speedMat;
  t1=nanmean(nanmean(tmpSM(3:6,:,:),3),1);
  e1 = nanstd(nanmean(tmpSM(3:6,:,:),3));%/sqrt(size(speedMatRaw,3));
t2=nanmean(nanmean(tmpSM(7:10,:,:),3),1);
e2 = nanstd(nanmean(tmpSM(7:10,:,:),3));%/sqrt(size(speedMatRaw,3));
figure('Position',[440   378   881   420],'Renderer','Painters')
subplot(2,2,1)
%errorbar(1:2:400,t1,e1)
boundedline(1:2:400,t1,e1,'cmap',[0 0 0],'alpha')
boundedline(1:2:400,t2,e2,'cmap',get_color(0.5,100),'alpha')
xticks([0:80:400])
for ii=1:4
    xline(ii*80,'k--')
end
% plot(t1,'k')
% hold on
% plot(t2,'Color',get_color(0.8,100))
% plot(t3,'k--')
subplot(2,2,3) 
t1=squeeze(nanmean(tmpSM(3:6,:,:)));
t2=squeeze(nanmean(tmpSM(7:10,:,:)));
pvals = zeros(1,200);
for ii=1:200
    pvals(ii)=signrank(t1(ii,:),t2(ii,:));
end
imagesc(repmat(pvals>=0.05/200,10,1))
axis image
colormap gray

subplot(2,2,2)
t1=squeeze(nanmean(nanmean(tmpSM(3:6,:,:),1),2));
t2=squeeze(nanmean(nanmean(tmpSM(7:10,:,:),1),2));
%plotSpread({t1-t2})
violin(t1-t2)
xlabel(sprintf('%.3f',signrank(t1,t2)))
saveas(gcf,fullfile('/Volumes/Samsung_T5/attialex/images','running05std.pdf'));
