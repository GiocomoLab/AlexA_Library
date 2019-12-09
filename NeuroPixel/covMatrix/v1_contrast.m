%% select datasets
sn = dir(fullfile(data_dir,'AA*_contr*.mat'));
for ii=1:size(sn,1)
    filenames{ii}=fullfile(sn(ii).folder,sn(ii).name);
end

%%  load datasets
contrast_levels = [100 50 20 10 5 2 0];
trials = 1:250;
region = 'VISp';
binsize=2;
contrast_stability = nan(numel(filenames),numel(contrast_levels));
contrast_stability_std = contrast_stability;
n_cells = zeros(1,numel(filenames));
for iF = 1:numel(filenames)
    data = load(filenames{iF});
    
    spatialMap=[];
    dwell_time=[];
    edges=[0:binsize:400];
    edges(1)=-.01;
    data.posx(data.posx<0)=0;
    data.posx(data.posx>=400)=399.00;
    for iT=1:length(trials)
        idxVR=data.trial==trials(iT);
        t_time=data.post(idxVR);
        start=min(t_time);
        stop=max(t_time);
        idxNP=data.sp.st<stop & data.sp.st>=start;
        [spM, dT]=getSpikeMatPosition(data.sp.st(idxNP),data.sp.clu(idxNP),data.posx(idxVR),data.post(idxVR),'edges',edges,'max_clust',max(data.sp.clu)+1);
        spatialMap=cat(3,spatialMap,spM);
        dwell_time=cat(1,dwell_time,dT);
    end
    %cellIDX=find(sp.cgs>=1);
    reg = strcmp(data.anatomy.cluster_parent,region);
    if iscolumn(reg)
        reg = reg';
    end
    spatialMap=spatialMap(data.sp.cids+1,:,:);
    spatialMap=spatialMap(data.sp.cgs==2 & reg,:,:);
    %spatialMap = spatialMap(this_stab>0.0 & this_MEC,:,:);
    %normalize by dwell time in each bin
    dt=dwell_time';
    dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
    for ii=1:size(spatialMap,1)
        spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
    end
    spatialMap(isnan(spatialMap))=0;
    % do spatial smoothing
    % filt = gausswin(11);
    % filt = filt/sum(filt);
    % filt = reshape(filt,[1, numel(filt),1]);
    smoothSigma = 4/binsize;
    smoothWindow = floor(smoothSigma*5/2)*2+1;
    gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
    filt = reshape(gauss_filter,[1, numel(gauss_filter),1]);
    sPF = repmat(spatialMap,[1,3,1]);
    sPF=convn(sPF,filt,'same');
    iidx = (size(spatialMap,2)+1):(2*size(spatialMap,2));
    sPF = sPF(:,iidx,:);
    % get template trials
    
    spatialMap = sPF;
    
    for iC=1:numel(contrast_levels)
        idx = data.trial_contrast == contrast_levels(iC) & data.trial_gain ==1;
        tmpMap = spatialMap(:,:,idx);
        stability = zeros(1,size(tmpMap,1));
        stabidx = find(triu(true(8),1));

for ii=1:numel(stability)
    tmp = squeeze(tmpMap(ii,:,:));
    tmp = corr(tmp);
    stability(ii)=mean(tmp(stabidx));
end
    contrast_stability(iF,iC)=nanmean(stability);
    contrast_stability_std(iF,iC)=nanstd(stability);
    end
    n_cells(iF)=size(spatialMap,1);
    
end

%%
set(0,'DefaultFigureRenderer','painters')


set(0,'DefaultFigureRenderer','painters')

%plot(contrast_stability','.-','MarkerSize',12)
errorbar(repmat(1:7,[4,1])',contrast_stability',contrast_stability_std')
set(gca,'XTick',[1:7],'XTickLabel',contrast_levels)
ylabel('Stability')
xlabel('Contrast Level [%]')
leg = {}
for ii=1:4
    leg{ii}=sprintf('Mouse %d, %d cells',ii,n_cells(ii));
end
legend(leg)
saveas(gcf,'C:/temp/contrast_stability_std.png')
saveas(gcf,'C:/temp/contrast_stability_std.pdf')

%%
set(0,'DefaultFigureRenderer','painters')

figure
%plot(contrast_stability','.-','MarkerSize',12)]
contrast_stability_sem =bsxfun(@rdivide, contrast_stability_std,sqrt(n_cells)');
errorbar(repmat(1:7,[4,1])',contrast_stability',contrast_stability_sem')
set(gca,'XTick',[1:7],'XTickLabel',contrast_levels)
ylabel('Stability')
xlabel('Contrast Level [%]')
leg = {}
for ii=1:4
    leg{ii}=sprintf('Mouse %d, %d cells',ii,n_cells(ii));
end
legend(leg)
saveas(gcf,'C:/temp/contrast_stability_sem.png')
saveas(gcf,'C:/temp/contrast_stability_sem.pdf')
    