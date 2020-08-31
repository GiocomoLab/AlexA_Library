matfiles = dir('/Volumes/Samsung_T5/attialex/tbtxcorr_decoder_time/*.mat');
pca_data = load('/Users/attialex/code/campbell_attinger/fig2_gain_response_types/pca_scores.mat');

errorMat = [];
distanceMat = [];
error_time_mat = [];
distance_time_mat = [];
take_idx_time = -100:100;
check_idx=1:100;
cluster_group = [];
figure
for iF=1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    cluster_this = pca_data.MEC_CLUSTERS(strcmp(pca_data.MEC_NAMES,sn));
    if isempty(cluster_this)
        continue
    end
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    cluster_group(end+1)=cluster_this;
    errorMat = cat(3,errorMat,squeeze(data_out.scoreMat(1,:,:)));
    distanceMat = cat(3,distanceMat,squeeze(data_out.scoreMat(2,:,:)));
    t1=data_out.time_error;
    t2=data_out.time_distance;
    
    remap_ons = strfind(abs(t1)>50,[zeros(1,60) ones(1,30)])+59;
    remap_trial = data_out.trials(remap_ons);
    if isempty(remap_trial)
        continue
    end
    remap_trial = remap_trial-(data_out.trials(1))+1;
    remap_ons=remap_ons(ismember(remap_trial,[7:10]));
    for iO=1:min(numel(remap_ons),30)
        
        
        error_sum=nnz(abs(t1(remap_ons(iO)+check_idx))>50);
        if error_sum>=75
            
            
            
            error_time_mat = cat(1,error_time_mat,t1(remap_ons(iO)+take_idx_time));
            distance_time_mat = cat(1,distance_time_mat,t2(remap_ons(iO)+take_idx_time));
        end
    end
    x_vec = 1:2:399;
    subplot(2,1,1)
    hold on
    for iT=1:16
        tmp = squeeze(data_out.scoreMat(1,iT,:));
        
        plot(x_vec+(iT-7)*400,tmp,'b')
        
    end
    axis tight
    ylim([-40 40])
    
    subplot(2,1,2)
    imagesc(squeeze(nanmean(data_out.corrMat)),[0 .7])
   % pause
    clf
end
%%
figure('Renderer','Painters')

subplot(2,1,1)
mm=mean((error_time_mat));
ee = std((error_time_mat))/sqrt(size(error_time_mat,1));
boundedline(take_idx_time/50,mm,ee)
grid on
xlabel('time from remap/failure')
ylabel('decoder error')
subplot(2,1,2)
tmp = distance_time_mat;%-mean(triggered_distance(:,10:20),2);
mm=nanmean(tmp);
ee = nanstd((distance_time_mat))/sqrt(size(distance_time_mat,1));
boundedline(take_idx_time/50,mm,ee,'alpha')
hold on

xlabel('time from remap/failure')
ylabel('similarity to baseline')
grid on
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_decoder_error_time.pdf'))
%%
errorMat_cluster = errorMat(:,:,cluster_group<3);
distanceMat_cluster = distanceMat(:,:,cluster_group<3);
dd=squeeze(mean(distanceMat_cluster(7,1:50,:)))-squeeze(mean(distanceMat_cluster(7,151:200,:)));
[~,sidx]=sort(dd);
slices = 2;
N=floor(size(errorMat_cluster,3)/slices);
figure('Renderer','Painters')
iT=7;
for iS=1:slices
    idx = (iS-1)*N+1:iS*N;
    idx=sidx(idx);
    subplot(2,1,1)
    hold on
    plot(squeeze(median(errorMat_cluster(iT,:,idx),3)))
    subplot(2,1,2)
    hold on
    plot(squeeze(median(distanceMat_cluster(iT,:,idx),3)));
end

%%
e_mat_cluster  = errorMat(:,:,cluster_group ~=3);
d_mat_cluster = distanceMat(:,:,cluster_group ~=3);
r = squeeze(mean(d_mat_cluster(7,1:50,:),2))-squeeze(mean(d_mat_cluster(7,150:200,:),2));
[~,sidx]=sort(r);
slices = 2;
N=floor(size(e_mat_cluster,3)/slices);
figure
for iS=1:slices
    idx = (iS-1)*N+1:iS*N;
    idx=sidx(idx);
    subplot(2,1,1)
    hold on
    plot(squeeze(median(e_mat_cluster(7,:,idx),3)))
    subplot(2,1,2)
    hold on
    plot(squeeze(median(d_mat_cluster(7,:,idx),3)))
end


%%
figure
subplot(2,1,1)
plot(-398:2:0,squeeze(median(errorMat(6,:,:),3)))
hold on
plot(2:2:400,squeeze(median(errorMat(7,:,:),3)))

subplot(2,1,2)
plot(-398:2:0,squeeze(median(distanceMat(6,:,:),3)))
hold on
plot(2:2:400,squeeze(median(distanceMat(7,:,:),3)))
%%
cmap = cbrewer('qual','Set2',3);
figure
x_vec = 1:2:399;
for iG=1:3
    cluster_idx = cluster_group==iG;
    subplot(2,1,1)
    hold on
    for iT=1:16
        plot(x_vec+(iT-7)*400,squeeze(median(errorMat(iT,:,cluster_idx),3)),'Color',cmap(iG,:))
        
    end
    subplot(2,1,2)
    hold on
    for iT=1:16
        plot(x_vec+(iT-7)*400,squeeze(median(distanceMat(iT,:,cluster_idx),3)),'Color',cmap(iG,:))
        
    end
end
%% trigger error on remap, 1 trigger per selected trial
triggered_error=[];
triggered_distance =[];
random_distance = [];
trial = 7;
take_idx = -20:20;
for ii=1:120
    t1=squeeze(errorMat(trial,:,ii));
    t2=squeeze(distanceMat(trial,:,ii));
    
    remap_ons = strfind(abs(t1)>50,[zeros(1,20) ones(1,20)])+19;
    if numel(remap_ons)>=1 && ismember(remap_ons(1),[21:180])
        rand_ons = randi([21,180],1);
        triggered_error = cat(1,triggered_error,t1(remap_ons(1)+take_idx));
        triggered_distance = cat(1,triggered_distance,t2(remap_ons(1)+take_idx));
        random_distance = cat(1,random_distance,t2(rand_ons+take_idx));
    end
    
end
figure('Renderer','Painters')
subplot(2,1,1)
mm=nanmean((triggered_error));
ee = std(abs(triggered_error))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee)

mm=mean(random_error);
ee = std((random_error))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee,'alpha','cmap',[.5 .5 .5])


grid on
xlabel('distance from remap/failure')
ylabel('error')
subplot(2,1,2)
tmp = triggered_distance;%-mean(triggered_distance(:,10:20),2);
mm=nanmean(tmp);
ee = std((triggered_distance))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee,'alpha')
hold on
mm=mean(random_distance);
ee = std((random_distance))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee,'alpha','cmap',[.5 .5 .5])
xlabel('distance from remap/failure')
ylabel('similarity to baseline')
grid on
%saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_decoder_error.pdf'))
%%
avg_dist = mean(distanceMat(7,1:50,:),2)-mean(distanceMat(7,151:end,:),2);

[~,sid]=sort(squeeze(avg_dist));

n_slices =4;
chunkSize = round(numel(avg_dist)/n_slices);
figure
for iSlice = 1:n_slices
    i_start = (iSlice-1)*chunkSize+1;
    i_stop = min(numel(avg_dist),iSlice*chunkSize);
    idx = sid(i_start:i_stop);
    subplot(2,1,1)
    hold on
    plot(squeeze(mean(distanceMat(7,:,idx),3)))
    subplot(2,1,2)
    hold on
    plot(squeeze(median(errorMat(7,:,idx),3)))
    
end



%%
%errorMatFlat=zeros(size(errorMat,1)*size(errorMat,2),size(errorMat,3));
%distanceMatFlat = errorMatFlat;
errorMatFlat = [];
distanceMatFlat = [];
for ii=7:10
    errorMatFlat = cat(1,errorMatFlat,squeeze(errorMat(ii,:,:)));
    distanceMatFlat = cat(1,distanceMatFlat,squeeze(distanceMat(ii,:,:)));
    
end
%%
triggered_error = [];
triggered_distance = [];
randomdistance=[];
random_error=[];
for ii=1:120
    t1=squeeze(errorMatFlat(:,ii))';
    t2=squeeze(distanceMatFlat(:,ii))';
    remap_ons = strfind(abs(t1)>50,[zeros(1,15),1 1 1 1 1 1])+14;
    for iO=1:min(numel(remap_ons),30)
        
        if ismember(remap_ons(iO),[21:(numel(t1)-20)])
            
            rand_ons = randi([21,(numel(t1)-20)],1);
            triggered_error = cat(1,triggered_error,t1(remap_ons(iO)+take_idx));
            triggered_distance = cat(1,triggered_distance,t2(remap_ons(iO)+take_idx));
            random_distance = cat(1,random_distance,t2(rand_ons+take_idx));
            random_error = cat(1,random_error,t1(rand_ons+take_idx));
        end
    end
end


figure('Renderer','Painters')
subplot(2,1,1)
mm=median((triggered_error*-1));
mmr=median(random_error*-1);
ee = mad((triggered_error))/sqrt(size(triggered_error,1));
eer = mad((random_error))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee)
boundedline(take_idx*2,mmr,eer,'alpha','cmap',[.5 .5 .5])

grid on
xlabel('distance from remap event')
ylabel('decoding error')
box off
subplot(2,1,2)
tmp = triggered_distance;%-mean(triggered_distance(:,10:20),2);
mm=mean(tmp);
ee = std((triggered_distance))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee,'alpha')
hold on
mm=mean(random_distance);
ee = std((random_distance))/sqrt(size(triggered_error,1));
boundedline(take_idx*2,mm,ee,'alpha','cmap',[.5 .5 .5])

grid on
xlabel('distance from remap event')
ylabel('cosine similarity')
box off
saveas(gcf,fullfile('/Users/attialex/Dropbox/temporary_images','fig6_decode_cos_sim.pdf'))
%%
figure
for ii=1:130
    %subplot(2,1,1)
    plot(-398:2:0,squeeze((errorMat(7,:,ii))),'b')
    hold on
    x_vec = 2:2:400;
    plot(x_vec,squeeze((errorMat(8,:,ii))),'b')
    ylim([-40 40])
    %     t1=squeeze(errorMat(7,:,ii));
    %     remap_ons = strfind(abs(t1)>50,[0 0 0 1 1 1])+3;
    %     if numel(remap_ons)>=1
    %         plot(x_vec(remap_ons(1)),t1(remap_ons(1)),'ro')
    %     end
    
    
    yyaxis right
    plot(-398:2:0,squeeze((distanceMat(7,:,ii))),'r')
    hold on
    plot(2:2:400,squeeze((distanceMat(8,:,ii))),'r')
    
    pause
    clf
end
