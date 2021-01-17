MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_random = MM_R-mean(MM_R(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
mm_resp = mean(MM_ms(:,opt.time_bins>0.1 & opt.time_bins<0.6),2);
mm_resp= mm_resp./FR;

valid_idx = FR>1;

%%
[~,sid]=sort(mm_resp(valid_idx));
MM_valid = MM_ms(valid_idx,:);
MM_random_valid = MM_random(valid_idx,:);
figure
imagesc(MM_random_valid(sid,:),[-1 1])
%%


%%
specgram_all=[];
for iP=1:numel(MM_THETA)
    specgram_all=cat(3,specgram_all,MM_THETA{iP}.mean_specgram);
end
specgram_valid = specgram_all(:,:,valid_idx);
theta_power = mean(mean(specgram_valid(s>6 & s<8,:,:),1),2);
[~,sid_theta]=sort(theta_power);
figure
imagesc(MM_valid(sid_theta,:),[-1 5])
%%
theta_power = mean(mean(specgram_valid(s>6 & s<8,:,:),1),2);
RANK=[];
[uS]=unique(SID(FR>1));
for iS=1:numel(uS)
    idx = SID(FR>1)==uS(iS);
    resp = theta_power(idx);
    ranking=1:numel(resp);
    [~,sid]=sort(resp);
    ranking(sid)=ranking/numel(resp);
    RANK = cat(1,RANK,ranking');
end

%%
%%
freq_idx = s>6 & s<8;
theta_power = mean(mean(specgram_valid(freq_idx,:,:),1),2);
theta_powerN = theta_power./mean(mean(specgram_valid(~freq_idx,:,:),1),2);
RANK=[];
[uS]=unique(SID(FR>1));
figure
hold on
data=[];
for iS=1:numel(uS)
    idx = SID(FR>1)==uS(iS);
    MM_this = MM_valid(idx,:);
    resp = theta_powerN(idx);
    ranking=1:numel(resp);
    [~,sid]=sort(resp);
    ranking(sid)=ranking/numel(resp);
    plot(mean(MM_valid(ranking<.1,:)))
    data = cat(1,data,mean(MM_valid(ranking>.9,:)));
end
figure
imagesc(data)
%%
figure
scatter(theta_powerN,templateDuration(FR>1),25,mm_resp(FR>1),'o','filled')
xlabel('theta power')
ylabel('template duration')
set(gca,'CLim',[-1 1])
colormap(brewermap(20,'*RdBu'))
%%
figure
scatter(theta_powerN,THETA_POWER(FR>1,1)-THETA_POWER(FR>1,2),'.')
xlabel('normalized theta power')
ylabel('template duration')

%%
quantiles = 10;
nSamples = size(MM_valid,1);
%cmap = cbrewer('div','RdBu');
cmap = brewermap(21,'RdBu');
cmap=flipud(cmap);
[~,sid_theta]=sort(theta_powerN);
xl =[-1 2.5];
for iC=[1 5 quantiles]%nChunks
    figure
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    sub_idx = round((iC-1)/quantiles*nSamples + 1):round(iC/quantiles*nSamples);
    IDX = sid_theta(sub_idx);
    %IDX=RANK>=(iC-1)/quantiles & RANK<iC/quantiles;
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    subplot(2,1,1)
    plotAVGSEM(MM_valid(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
    plotAVGSEM(MM_random_valid(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])
    
    xlim(xl)
    
    subplot(2,1,2)
    imagesc(opt.time_bins,1:nnz(IDX),MM_valid(IDX,:),[-1 5])
    xlim(xl)
    %colormap(cmap);
    
end
%%
%%
quantiles = 10;
nSamples = size(MM_valid,1);
%cmap = cbrewer('div','RdBu');
cmap = brewermap(21,'RdBu');
cmap=flipud(cmap);
xl =[-1 2.5];
Waveform_valid = templateWaveform(FR>1,:);
[~,sid_duration]=sort(templateDuration(FR>1));

for iC=[1 5 quantiles]%nChunks
    figure
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    sub_idx = round((iC-1)/quantiles*nSamples + 1):round(iC/quantiles*nSamples);
    IDX = sid_duration(sub_idx);
    %IDX=RANK<.1
    %IDX=RANK>=(iC-1)/quantiles & RANK<iC/quantiles;
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    subplot(2,1,1)
    plotAVGSEM(smoothdata(MM_valid(IDX,:),2,'gaussian',7)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
    plotAVGSEM(smoothdata(MM_random_valid(IDX,:),2,'gaussian',7)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])
    
    xlim(xl)
    
    subplot(2,1,2)
    boundedline(1:82,mean(Waveform_valid(IDX,:)),std(Waveform_valid(IDX,:)))
    figure
    imagesc(opt.time_vecs,1:numel(IDX),smoothdata(MM_valid(IDX,:),2,'gaussian',11),[-1 3])

end


%%
figure
waterfall(s(s>1),time_vec,mean(specgram_valid(s>1,:,:),3)')
ylim([-5 5])
%%
figure
waterfall(s(s>1),time_vec,mean(specgram_valid(s>1,:,sid(end-(1:100))),3)')
ylim([-5 5])
%%
figure
[~,sid]=sort(templateDuration(FR>1));
%[~,sid]=sort(theta_powerN);
subplot(1,2,1)
imagesc(Waveform_valid(sid,:),[-.4 .4])
xlim([30,70])
title('template waveform')
subplot(1,2,2)
imagesc(opt.time_vecs,1:numel(sid),smoothdata(MM_valid(sid,:),2,'gaussian',11),[-1 3])
title('MM response')
xlim([-1 3])