%matfiles =
%dir('/Volumes/Samsung_T5/attialex/tbtxcorr_decoder_0.5_fitcoec_svm/*.mat');
%% with speed cutoff

%% without speed cutoff
matfiles = dir('/Volumes/Samsung_T5/attialex/tbtxcorr_decoder_0.5_noSpeedCutoff/*.mat');

errorMat = [];
distanceMat = [];
error_time_mat = [];
error_time_mat_2 = [];
ALLOnsets=[];
remap_time = [];
error_classifier = [];
distance_time_mat = [];
take_idx_time = -450:450;
take_idx_space = -20:20;
check_idx=1:100;
x_vec = 1:2:399;
ops = load_default_opt;
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.xbinedges = ops.edges;
ops.xbincent = .5*ops.edges(1:end-1)+.5*ops.edges(2:end);

for iF=1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    errorMat = cat(3,errorMat,squeeze(data_out.scoreMat(1,:,:)));
    distanceMat = cat(3,distanceMat,squeeze(data_out.scoreMat(2,:,:)));
    t1=data_out.time_error*-1;
    t2=data_out.time_distance;
    
    remap_ons = strfind(abs(t1)>50,[zeros(1,150) 1])+151;
    remap_trial = data_out.trials(remap_ons);
    
    remap_trial = remap_trial-(data_out.trials(1))+1;
    remap_ons=remap_ons(ismember(remap_trial,[7]));
    if isempty(remap_ons)
        continue
    end
    iO=1;
    trials = data_out.trials-data_out.trials(1)+1;
    
    gain_trial_onset = strfind(trials'==7,[0 1]);
    remap_time = cat(1,remap_time,remap_ons(iO)-gain_trial_onset);
    
    error_time_mat = cat(1,error_time_mat,t1(remap_ons(iO)+take_idx_time));
    
    
    remap_trial = data_out.trials;
    
    remap_trial = remap_trial-(data_out.trials(1))+1;
    remap_ons = strfind((remap_trial)'==7,[zeros(1,20) ones(1,20)])+19;
    [~,~,posbin] = histcounts(data_out.posx,ops.xbinedges);
    
    for iO=1:min(numel(remap_ons),1)
        error_time_mat_2 = cat(1,error_time_mat_2,t1(remap_ons(iO)+take_idx_time));
    end
end


%% plot trials triggered on remap
val = mad((error_time_mat(:,take_idx_time<-2)),1,2);
[~,sid]=sort(val);
%[~,sid]=sort(remap_time);
figure('Renderer','Painters','Position',[440   442   693   356])
subplot(2,2,1)
imagesc(take_idx_time/50,1:size(error_time_mat,1),error_time_mat(sid,:),[-10 10])
N=size(error_time_mat,1);
%N=70;
subplot(2,2,2)
hold on

tmp = squeeze(nanmedian(error_time_mat(sid(1:N),:)));
tmpMAD = mad(error_time_mat(sid(1:N),:),1);
%plot(x_vec+(iT-7)*400,tmp,'b')
%errorbar(x_vec+(iT-7)*400,tmp,tmpMAD,'b')
boundedline(take_idx_time/50,tmp,tmpMAD,'b','alpha')
xlim([-8 0.1])
ylim([-10 30])
%subplot(2,2,4)
xlabel('time to remap')
ylabel('classifier error')
box off
subplot(2,2,3)
plotSpread(-1*remap_time(sid(1:N))/50,'xyOri','flipped','spreadWidth',.1,'spreadFcn',{'xp',.1})
hold on
boxplot(remap_time*-1/50,'PlotStyle','compact','Orientation','horizontal')
xlim([-8 0.1])
box off


N=[1 10 24];
ind = 0;
%N=1:40;
cmap = cbrewer('qual','Set2',numel(N));

for ii=N
    ind = ind+1;
    subplot(2,2,4)
    hold on
    
    tmp = error_time_mat(sid(ii),:);
    remap_ons = strfind(abs(tmp)>=50& take_idx_time>0,[0 0 0 1])+1;
    tmp(remap_ons(1):end)=nan;
    tmp2 = error_time_mat(sid(ii),:);
    tmp2(1:remap_ons(1)-2)=nan;
    
    plot(take_idx_time/50,tmp,'Color',cmap(ind,:),'LineWidth',2)
    xline(-remap_time(sid(ii))/50,'Color',cmap(ind,:),'LineWidth',2);
    plot(take_idx_time/50,tmp2,'--','Color',cmap(ind,:))
    remap_time(sid(ii))
    title(ii)
    %pause
    %clf
end
ylim([-80 80])
xlim([-7 7])
box off
xlabel('time to remap [s]')
ylabel('classifier error')
saveas(gcf,'/Users/attialex/Dropbox/temporary_images/fig6_remap_time.pdf')



%% triggered on gain change onset
val = mean(abs(error_time_mat_2(:,take_idx_time<0)),2);
[~,sid]=sort(val);
figure
subplot(2,1,1)
N=size(error_time_mat_2,1);
error_mat_corrected = nan(N,numel(take_idx_time));
tmp_times= [];
for ii=1:N
    tmp = error_time_mat_2(sid(ii),:);
    remap_ons = strfind(abs(tmp)>=50& take_idx_time>0,[0 0 0 1])+3;
    error_preceeding = tmp(take_idx_time/50>-5 & take_idx_time/50<0);
    if nnz(abs(error_preceeding)>50)
        %continue
    end
    tmp_times(end+1)=remap_ons(1);
    tmp(remap_ons(1):end)=nan;
    hold on
    plot(take_idx_time/50,tmp,'Color',[.8 .8 .8])
    error_mat_corrected(ii,:)=tmp;
end
hold on
plot(take_idx_time/50,nanmedian(error_mat_corrected),'k','LineWidth',2)
xlabel('time from gain onset')
%xlim([-1.5 6])
xlim([-2 8.5])

box off
ylim([-40 40])
xlabel('time from gain change onset')
ylabel('classifier error')
subplot(2,1,2)
hold on
tmp_times = tmp_times/50+take_idx_time(1)/50;
plotSpread(tmp_times','xyOri','flipped','spreadWidth',.1,'spreadFcn',{'xp',.1})
boxplot(tmp_times,'PlotStyle','compact','Orientation','horizontal')
xlim([-2 8.5])
xlabel('time from gain change onset')

box off

%saveas(gcf,'/Users/attialex/Dropbox/temporary_images/fig6_remap_trial_onset.pdf')
%%
%% plot median and median absolute deviation trial onset
figure
subplot(1,2,1)
hold on
for iT=1:16
    tmp = squeeze(nanmedian(errorMat(iT,:,:),3));
    tmpMAD = squeeze(mad(errorMat(iT,:,:),1,3));
    %plot(x_vec+(iT-7)*400,tmp,'b')
    %errorbar(x_vec+(iT-7)*400,tmp,tmpMAD,'b')
    boundedline(x_vec+(iT-7)*400,tmp,tmpMAD)
    
end

subplot(1,2,2)
hold on

tmp = squeeze(nanmedian(error_time_mat_2));
tmpMAD = mad(error_time_mat_2,1);
%plot(x_vec+(iT-7)*400,tmp,'b')
%errorbar(x_vec+(iT-7)*400,tmp,tmpMAD,'b')
boundedline(take_idx_time/50,tmp,tmpMAD)


%% error time mat 2 for session with small MAD preceeding
val = mean(abs(error_time_mat_2(:,take_idx_time<0)),2);
[~,sid]=sort(val);
figure
subplot(1,3,1)
imagesc(take_idx_time/50,1:size(error_time_mat_2,1),error_classifier(sid,:),[-10 10])

subplot(1,3,2)
hold on
N=size(error_time_mat_2,1);
% N=40;
tmp = squeeze(nanmedian(error_time_mat_2(sid(1:N),:)));
tmpMAD = mad(error_time_mat_2(sid(1:N),:),1);
%tmp = squeeze(nanmean(error_time_mat_2(sid(1:N),:)));
%tmpMAD = squeeze(nanstd(error_time_mat_2(sid(1:N),:)))/sqrt(N);
%plot(x_vec+(iT-7)*400,tmp,'b')
%errorbar(x_vec+(iT-7)*400,tmp,tmpMAD,'b')
boundedline(take_idx_time/50,tmp,tmpMAD,'b','alpha')

subplot(1,3,3)
hold on


tmp = squeeze(nanmedian(error_classifier(sid(1:N),:)));
tmpMAD = mad(error_classifier(sid(1:N),:),1);
%tmp = squeeze(nanmean(error_time_mat_2(sid(1:N),:)));
%tmpMAD = squeeze(nanstd(error_time_mat_2(sid(1:N),:)))/sqrt(N);
%plot(x_vec+(iT-7)*400,tmp,'b')
%errorbar(x_vec+(iT-7)*400,tmp,tmpMAD,'b')
boundedline(take_idx_time/50,tmp,tmpMAD,'b','alpha')
figure
N=10;
cmap = cbrewer('div','Spectral',N);
for ii=1:N
    subplot(1,2,1)
    hold on
    plot(take_idx_time/50,error_classifier(sid(ii),:)','Color',cmap(ii,:))
    title('Classifier')
    subplot(1,2,2)
    hold on
    plot(take_idx_time/50,error_time_mat_2(sid(ii),:)')
    title('Distance based')
end

