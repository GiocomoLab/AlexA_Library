%load('F:\G5\1207_mismatch_1\1207_mismatch_1.mat')
PLOT=true;
speed_t=0.05;
% figure('Name',filenames{iF});; plot(speed)
%
if size(mismatch_trigger,1) ~=1
    mismatch_trigger=mismatch_trigger';
end
speed=true_speed';

all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
run_periods=smooth(speed,25)>speed_t;
run_window=-30:30;
possibles=strfind(run_periods',ones(1,length(run_window)))+floor(.5*length(run_window));


mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
possibles=randsample(possibles,500);
%% MM
[spike_mat,win,adata]=extract_triggered_spikes(sp,post(mm_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
[spike_mat_random,~,adata_random]=extract_triggered_spikes(sp,post(possibles),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
[spike_mat_all,~,adata_all]=extract_triggered_spikes(sp,post(all_mm_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
%% RO
ron_trigs=strfind(run_periods',[zeros(1,100),ones(1,100)])+100;
roff_trigs = strfind(run_periods',[ones(1,100),zeros(1,100)])+100;

valid_idx=true(size(roff_trigs));
for ii=1:length(roff_trigs)
    tr=roff_trigs(ii);
    rng=tr+(-50:1:10);
    closeMMTrigs=ismember(all_mm_trigs,rng);
    if nnz(closeMMTrigs>0)
        valid_idx(ii)=0;
    end
end
roff_trigs=roff_trigs(valid_idx);
%%
[spike_matRON,win,adataRON]=extract_triggered_spikes(sp,post(ron_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
[spike_matROFF,win,adataROFF]=extract_triggered_spikes(sp,post(roff_trigs),'win',[-4 4],'aux',[post'; [speed; mismatch_trigger]],'aux_win',[-200 200]);

% figure('Name',filenames{iF});
% subplot(1,2,1)
% imagesc(squeeze(adataROFF(1,:,:)));
% subplot(1,2,2)
% imagesc(squeeze(adataROFF(2,:,:)));
%% delete the rows containing no data
spike_mat=spike_mat(sp.cids+1,:,:);
spike_mat_random=spike_mat_random(sp.cids+1,:,:);
spike_mat_all = spike_mat_all(sp.cids+1,:,:);

spike_matRON = spike_matRON(sp.cids+1,:,:);
spike_matROFF = spike_matROFF(sp.cids+1,:,:);


%% turn spike_mat into firing rate
kernel=reshape(gausswin(401),1,1,[]);
binsize=0.001;
kernel=kernel/sum(kernel(:))/binsize;
%kernel=gpuArray(kernel);
mm_rate=convn(spike_mat,kernel,'same');
%%
zn_random = convn(spike_mat_random,kernel,'same');
rON_rate = convn(spike_matRON,kernel,'same');
rOFF_rate  = convn(spike_matROFF,kernel,'same');
%% separate into matrix for plotting and sorting, plus use only clusters with certain ID
cellIDX=sp.cgs>=2;
%cellIDX=true(1,size(zn,1));
mean_fr_sort=squeeze(mean(mm_rate(cellIDX,2:2:end,:),2));
mean_fr_plot=squeeze(mean(mm_rate(cellIDX,1:2:end,:),2));
mean_fr_rand=squeeze(mean(zn_random(cellIDX,:,:),2));
mean_fr_roff=squeeze(mean(rOFF_rate(cellIDX,:,:),2));

%% plot and sort cells
if PLOT
    [~,max_c]=max(mean_fr_sort,[],2);
    mm_resp=mean(mean_fr_sort(:,4100:4599),2)-mean(mean_fr_sort(:,3500:3999),2);
    [ss,sort_idx]=sort(mm_resp,'desc');
    %[ss,ii]=sort(max_c);
    
    
    [a,~]=max(mean_fr_plot,[],2);
    ff_plot=bsxfun(@rdivide,mean_fr_plot,a);
    ff_sort=bsxfun(@rdivide,mean_fr_sort,max(mean_fr_sort,[],2));
    
    %ff_rand=bsxfun(@rdivide,mean_fr_rand,max(mean_fr_rand,[],2));
    %ff_roff=bsxfun(@rdivide,mean_fr_roff,max(mean_fr_roff,[],2));
    
    %figure('Name',filenames{iF});('Name',filenames{iF})
    try
    figure('Name',filenames{iF});
    catch ME
        figure
    end
    subplot(3,1,1)
    imagesc(ff_sort(sort_idx,:))
    title('data used for sorting')
    subplot(3,1,2)
    imagesc(ff_plot(sort_idx,:))
    title('held out data')
    subplot(3,1,3)
    plot(-4:0.02:4,squeeze(mean(adata)))
    hold on
    plot(-4:0.02:4,squeeze(mean(adata_random)))
    
    plot(-4:0.02:4,squeeze(mean(adataROFF(1,:,:))),'k-')
    try
    figure('Name',filenames{iF});
    catch
        figure;
    end
    imagesc(ff_roff(sort_idx,:))
    title('Run Offset Response')
    
    %%
    rOffresp=mean(mean_fr_roff(:,4001:4500),2)-mean(mean_fr_roff(:,3251:3750),2);
    [a,b]=sort(rOffresp);
    try
    figure('Name',filenames{iF});
    catch
        figure;
    end
    imagesc(ff_roff(b,:))
    title('Run Offset Response')
    mean_fr_mm=squeeze(mean(mm_rate(cellIDX,:,:),2));
    mean_mm=mean(mean_fr_mm(:,4101:4599),2)-mean(mean_fr_mm(:,3500:3999),2);
    plot(mean_mm,rOffresp,'.')
    lims=[-5 5];
    xlim(lims);
    ylim(lims)
    axis image
    grid on
    ylabel('offset response')
    xlabel('mm response')
    %%
    nspikesR=sum(spike_mat_random(:,:,4001:5000),3);
    nspikes=sum(spike_mat(:,:,4001:5000),3);
    
    SIG=zeros(size(nspikesR,1),1);
    for ii=1:size(nspikesR,1)
        SIG(ii)=ranksum(nspikes(ii,:),nspikesR(ii,:));
    end
    tmp=SIG(cellIDX);
    try
    figure('Name',filenames{iF});
    catch
        figure;
    end
    imagesc(tmp(sort_idx)<0.05);
    
    %%
    resp = squeeze(mean(mm_rate(sp.cgs==2,:,:),2));
    respR = squeeze(mean(zn_random,2));
    roffResp = squeeze(mean(rOFF_rate,2));
    ronResp = squeeze(mean(rON_rate,2));
    
    params=struct();
    params.winIDX=-4000:4000;
    params.masterTime=params.winIDX/1000;
    params.xLim=[-1 3];
    paramsAdata=struct();
    paramsAdata.winIDX=-200:200;
    paramsAdata.masterTime=paramsAdata.winIDX/50;
    paramsAdata.xLim=[-1 3];
    try
    figure('Name',filenames{iF});
    catch
        figure
    end
    subplot(2,1,1)
    plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999)
    
    plotAVGSEM(respR',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',[1 0 0]);
    legend({'MM','Random'})
    
    subplot(2,1,2)
    plotAVGSEM(roffResp',gca,'parameters',params,'ms',true,'baseline',3500:3999)
    plotAVGSEM(ronResp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',[1 0 0])
    legend({'Run Off','Run On'})
    %%
    figure('Name',filenames{iF});
    subplot(211)
    plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999)
    plotAVGSEM(roffResp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',[.8 .8 0])
    plotAVGSEM(ronResp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',[1 0 0])
    legend({'MM','Run Off','Run On'});
    xlabel('time [s]')
    ylabel('Change in Firing Rate [Hz]')
    
    
    subplot(2,1,2)
    plotAVGSEM(squeeze(adata)',gca,'parameters',paramsAdata,'ms',false)
    hold on
    plotAVGSEM(squeeze(adataROFF(1,:,:))',gca,'parameters',paramsAdata,'ms',false,'col',[.8 .8 0])
    plotAVGSEM(squeeze(adataRON(1,:,:))',gca,'parameters',paramsAdata,'ms',false,'col',[1 0 0])
    plotAVGSEM(squeeze(adata_random)',gca,'parameters',paramsAdata,'ms',false,'col',[0 0 0])
    legend({'MM','Run Off','Run On','Random'});
    title('Run Traces')
    xlabel('time [s]')
    ylabel('Run speed')
    
    %%
    figure('Name',filenames{iF});
    shift = 0.76;
    subplot(2,1,1)
    plot(-4:0.02:4,squeeze(mean(adata)))
    hold on
    plot([-4:0.02:4]+shift,squeeze(mean(adataROFF(1,:,:))))
    
    subplot(2,1,2)
    plot(-4:0.001:4,mean(resp))
    hold on
    plot([-4:0.001:4]+shift,mean(roffResp))
    legend({'MM','Shifted rOff'})
    
    %%
    cell_list=find(cellIDX);
    cellID=cell_list(sort_idx(2));
    binned_array=squeeze(spike_mat(cellID,:,:));
    bins=win(1):0.001:win(2);
    
    % set the new rasters
    [tr,b] = find(binned_array);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy*1+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
    figure('Name',filenames{iF});;subplot(2,1,1);plot(rasterX,rasterY)
    hold on
    
    subplot(2,1,2)
    plot(-4:0.001:4,resp(cellID,:));
    hold on
    plot(-4:0.001:4,respR(cellID,:))
    %%
    % figure('Name',filenames{iF});
    % plot(squeeze(mean(adata,2)))
    % speed_diff=mean(adata(:,:,225:275),3)-mean(adata(:,:,150:200),3);
    % hold on
    % [a,b]=sort(speed_diff);
    % plot(squeeze(mean(adata(:,b(1:20),:),2)))
    % plot(squeeze(mean(adata(:,b(end-20:end),:),2)))
    
    runS=mean(adata_all(:,:,175:225),3);
    [~,~,cOff]=unique(runS);
    rank=cOff/max(cOff)';
    nSplit=4;
    frac=1/nSplit;
    kk=reshape(gausswin(401),1,[]);
    figure('Name',filenames{iF});
    col=summer(4);
    for ii=1:nSplit
        lims=frac*[ii-1 ii];
        tmpidx=rank>lims(1) & rank<=lims(2);
        
        resp = squeeze(mean(spike_mat_all(sp.cgs==2,tmpidx,:),2));
        resp = convn(resp,kk,'same');
        subplot(2,1,1)
        plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',col(ii,:))
        subplot(2,1,2)
        hold on
        plot(-4:0.02:4,squeeze(mean(adata_all(:,tmpidx,:),2)),'Color',col(ii,:))
    end
    
end