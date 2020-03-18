ops.BinWidth = 2;

ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.TimeBin = 0.02;
region = 'VISp';
data_path = 'F:/NP_DATA';
trgain = [ones(1,6), 0.8*ones(1,4), ones(1,6)];

%%
%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
matfiles = dir('F:\temp\tbtxcorr\*.mat');
allM = [];
allS = [];
allSpeed = [];
allSpeedZ= [];
SUBMATAdjusted = [];
SUBMAT = [];
baseline_shifts = [];
similarity_score = [];
similarity_score_MEC=[];
cntr = 0;
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    idx_region  = startsWith(data_out.region,region);
    idx = data_out.stability>.4 & startsWith(data_out.region,region)';
    if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.1
        
        cntr = cntr+1;
        tmp=squeeze(nanmean(data_out.corrMat(idx,:,:)));
        ff=tmp(7:16,1:6);
        similarity_score(end+1) = mean(ff(:));
        
        tmpS = squeeze(nanmean(data_out.shiftMat(idx,:,:)));
        allM=cat(3,allM,tmp);
        allS = cat(3,allS,tmpS);
        %baseline_shifts = cat(1,baseline_shifts,reshape(tmp_peak',1,[]));
        sn = [matfiles(iF).name(1:end-6) '.mat'];
        data = load(fullfile(data_path,sn),'posx','post','trial_gain','trial');
        [speed,speed_raw]=calcSpeed(data.posx,ops);
        nT=numel(data.trial_gain);
        speedMatAdjusted = nan(nT,ops.nBins);
        speedMat = nan(nT,ops.nBins);
        
        
        discrete_pos = discretize(data.posx,ops.edges);
        
        for iT=data_out.trials
            for iB=1:ops.nBins
                idx = data.trial==iT & discrete_pos == iB;
                if nnz(idx)>0
                    val = mean(speed(idx));
                    speedMatAdjusted(iT,iB)=val/data.trial_gain(iT);
                    speedMat(iT,iB)=val;
                end
            end
            if nnz(isnan(speedMatAdjusted(iT,:)))>0
            temp = speedMatAdjusted(iT,:);
            temp2 = speedMat(iT,:);
            nanx=isnan(temp);
            t=1:numel(temp);
            speedMatAdjusted(iT,isnan(temp)) = interp1(t(~nanx),temp(~nanx),t(nanx),'pchip');
            speedMat(iT,isnan(temp2)) = interp1(t(~nanx),temp2(~nanx),t(nanx),'pchip');
            end
        end
        subMatAdjusted=speedMatAdjusted(data_out.trials,:);
        SUBMATAdjusted = cat(3,SUBMATAdjusted,subMatAdjusted);
        SUBMAT = cat(3,SUBMAT,speedMat(data_out.trials,:));
        cc=corr(subMatAdjusted');
        cc=cc-diag(diag(cc));
        allSpeed = cat(3,allSpeed,cc);
        subMatZ = zscore(subMatAdjusted,1,2);
        dd=corr(subMatZ');
        dd=dd-diag(diag(dd));
        allSpeedZ = cat(3,allSpeedZ,dd);
    else
        nstab = nnz(idx);
        ntot = nnz(idx_region);
        sprintf('Kicked %s, stab: %d, tot: %d \n',matfiles(iF).name,nstab,ntot)
    end
end
%%
num_tr = size(allM,1);
tf = triu(true(num_tr),1);
Y = [];
for ii=1:size(allM,3);
    tmp = allM(:,:,ii);
    Y=cat(1,Y,tmp(tf)');
end
[coeff,score,~,~,expl] = pca(Y);
[~,sort_idx_gc] = sort(score(:,1),'descend');
corrmat_sort = allM(:,:,sort_idx_gc);

savefig{2} = figure('Position',[100 100 1000 800],'Renderer','Painters');
ha = tight_subplot(10,10);

savefig{2} = figure('Position',[100 100 1000 800],'Renderer','Painters');
ha2 = tight_subplot(10,10);

for i = 1:size(corrmat_sort,3)
    axes(ha(i)); hold on;
    imagesc(squeeze(corrmat_sort(:,:,i)));
    caxis([0 0.7]);
    %colorbar;
    %axis square;
    %title(sprintf('session %d',i));
    
    % patches indicating gain value
    for tr = 1:num_tr
        patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
        patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
    end
    axis image
    
    axes(ha2(i)); hold on;
    imagesc(squeeze(allSpeedZ(:,:,sort_idx_gc(i))));
    caxis([0 0.7]);
    %colorbar;
    %axis square;
    %title(sprintf('session %d',i));
    
    % patches indicating gain value
    for tr = 1:num_tr
        patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
        patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
            'EdgeColor',get_color(trgain(tr),100));
    end
    axis image

end

%%
subMat_normalized = SUBMATAdjusted;
for iS=1:size(SUBMATAdjusted,3)
    fact = max(SUBMATAdjusted(:,:,iS),[],'all');
    subMat_normalized(:,:,iS)=subMat_normalized(:,:,iS)/fact;
end

%%
num_tr = size(allM,1);
tf = triu(true(num_tr),1);
Y = [];
for ii=1:size(allM,3)
    tmp = allM(:,:,ii);
    Y=cat(1,Y,tmp(tf)');
end
[coeff,score,~,~,expl] = pca(Y);
[~,sort_idx_gc] = sort(score(:,1),'descend');
corrmat_sort = allM(:,:,sort_idx_gc);


figure
subplot(3,2,1)
plot(squeeze(mean(subMat_normalized(6,:,sort_idx_gc(1:20)),3)))
hold on
plot(squeeze(mean(subMat_normalized(7,:,sort_idx_gc(1:20)),3)))
plot(squeeze(mean(subMat_normalized(8,:,sort_idx_gc(1:20)),3)))
ylim([0 1])
legend('BL','G1','G2')
subplot(3,2,2)
plot(squeeze(mean(subMat_normalized(6,:,sort_idx_gc(end-19:end)),3)))
hold on
plot(squeeze(mean(subMat_normalized(7,:,sort_idx_gc(end-19:end)),3)))
plot(squeeze(mean(subMat_normalized(8,:,sort_idx_gc(end-19:end)),3)))
ylim([0 1])
legend('BL','G1','G2')

subplot(3,2,3)
imagesc(mean(allM(:,:,sort_idx_gc(1:10)),3),[.3 0.65])
axis image
title('Firing Rate Xcorr')
subplot(3,2,4)
imagesc(mean(allM(:,:,sort_idx_gc(end-9:end)),3),[.3 0.65])
axis image
title('Firing Rate Xcorr')

subplot(3,2,5)
imagesc(mean(allSpeed(:,:,sort_idx_gc(1:10)),3),[.3 0.65])
axis image
title('Run Speed Xcorr')

subplot(3,2,6)
imagesc(mean(allSpeed(:,:,sort_idx_gc(end-9:end)),3),[.3 0.65])
axis image
title('Run Speed Xcorr')
%%
averageShift = squeeze(mean(mean(allS(1:6,7:10,:),1),2));
averageSpeedPre = mean(mean(SUBMATAdjusted(1:6,:,:),1),2);
averageSpeedGC = mean(mean(SUBMATAdjusted(7:10,:,:),1),2);
figure
plot(squeeze(averageSpeedPre-averageSpeedGC),averageShift,'.')

averageShift = squeeze(mean(mean(allS(1:6,7:10,:),1),2));
averageSpeedPre = mean(mean(SUBMAT(1:6,:,:),1),2);
averageSpeedGC = mean(mean(SUBMAT(7:10,:,:),1),2);
figure
plot(squeeze(averageSpeedPre-averageSpeedGC),averageShift,'.')