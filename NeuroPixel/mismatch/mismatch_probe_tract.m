MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_random = MM_R-mean(MM_R(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
mm_resp = mean(MM_ms(:,opt.time_bins>0.1 & opt.time_bins<0.6),2);
mm_resp= mm_resp./FR;

valid_idx = FR>1 & ~isnan(POS3D(:,3));
[~,sid]=sort(mm_resp(valid_idx));
MM_valid = smoothdata(MM_ms(valid_idx,:)./FR(valid_idx),2,'gaussian',5);
MM_random_valid = smoothdata(MM_random(valid_idx,:)./FR(valid_idx),2,'gaussian',5);

%%
figure
h=scatter3((POS3D(valid_idx,1)),POS3D(valid_idx,2),POS3D(valid_idx,3),15,mm_resp(valid_idx),'filled');
xlabel('ml')
zlabel('dv')
axis equal
h.MarkerEdgeAlpha=0.5;
h.MarkerFaceAlpha = 0.5;
set(gca,'CLim',[-1 1])
colormap(brewermap(12,'*RdBu'))

%%
specgram_all=[];
for iP=1:numel(MM_THETA)
    specgram_all=cat(3,specgram_all,MM_THETA{iP}.mean_specgram);
end
freq_idx = s>6 & s<8;
theta_power = mean(mean(specgram_all(freq_idx,:,:),1),2);
theta_powerN = squeeze(theta_power./mean(mean(specgram_all(~freq_idx,:,:),1),2));

figure
h=scatter3((POS3D(valid_idx,1)),POS3D(valid_idx,2),POS3D(valid_idx,3),15,theta_powerN(valid_idx),'filled')
xlabel('ml')
zlabel('dv')
axis equal
set(gca,'CLim',[0 2])
colormap(brewermap(12,'*RdBu'))

%%
quantiles = 10;
nSamples = size(MM_valid,1);
%cmap = cbrewer('div','RdBu');
cmap = brewermap(21,'RdBu');
cmap=flipud(cmap);
[~,sid_z]=sort(POS3D(valid_idx,3),'descend');
xl =[-1 2.5];
for iC=[1 5 quantiles]%nChunks
    figure
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    sub_idx = round((iC-1)/quantiles*nSamples + 1):round(iC/quantiles*nSamples);
    IDX = sid_z(sub_idx);
    %IDX=RANK>=(iC-1)/quantiles & RANK<iC/quantiles;
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    subplot(2,1,1)
    plotAVGSEM(MM_valid(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0)
    plotAVGSEM(MM_random_valid(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',[.5 .5 .5])
    
    xlim(xl)
    
    subplot(2,1,2)
    imagesc(opt.time_bins,1:nnz(IDX),MM_valid(IDX,:),[-1 1])
    xlim(xl)
    colormap(brewermap(13,'*RdBu'))
    %colormap(cmap);
    
end
%%
figure
quantiles = 5;
nSamples = size(MM_valid,1);
%cmap = cbrewer('div','RdBu');
cmap = brewermap(quantiles,'RdYlBu');
[~,sid_z]=sort(POS3D(valid_idx,3),'descend');
xl =[-1 2.5];
for iC=1:quantiles%nChunks
    
    sub_idx=(iC-1)*chunksize+(1:chunksize);
    sub_idx = round((iC-1)/quantiles*nSamples + 1):round(iC/quantiles*nSamples);
    IDX = sid_z(sub_idx);
    %IDX=RANK>=(iC-1)/quantiles & RANK<iC/quantiles;
    params=struct();
    params.masterTime=opt.time_bins(1:end-1)*.5+opt.time_bins(2:end);
    params.xLim=[-2 3];
    hold on
    plotAVGSEM(MM_valid(IDX,:)',gca,'parameters',params,'ms',false,'baseline',opt.time_bins>=-.5 & opt.time_bins<0,'col',cmap(iC,:))
    
    xlim(xl)
    

    
end
%%
figure
imagesc(opt.time_vecs,1:nnz(valid_idx),MM_valid(sid_z,:),[-1 1])
colormap(brewermap(12,'*RdBu'))
xlim([-1 2.5])