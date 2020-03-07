ops.factors = -.55:0.01:.55;
ops.BinWidth =1;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.n_preceeding = 18;
ops.trials = 3:20;
ops.TimeBin = 0.02;
ops.idx = [10:ops.BinWidth:390]/ops.BinWidth;% in bins
fi = gausswin(11);
fi=fi'/sum(fi);
ops.filter = fi;
ops.bl_pre = 1:ops.n_preceeding;
ops.gain_trials = ops.n_preceeding+[1:4];
ops.similarity_trials = ops.n_preceeding+[-5:4];
ops.plotfig = false;
ops.maxLag = 20; % in cm
%%
gain = 0.8;
contrast = 100;
region = 'VISp';
[filenames,triggers] = getFilesCriteria(region,contrast,gain,'F:/NP_DATA');
%%
p=gcp('nocreate');
if isempty(p)
    parpool();
end
shift_data_path = 'Z:\giocomo\attialex\speed_filtered_smallBin';
%%
output = cell(numel(filenames),1);
parfor iF=1:numel(filenames)
    data = load(filenames{iF});
    ops_local = ops;
    data_out=cell(2,1);
    for iRep = 1:min(2,numel(triggers{iF}))
        preceeding_trials = triggers{iF}(iRep)+[-ops_local.n_preceeding:-1];
        if ~all(data.trial_gain(preceeding_trials)==1 & data.trial_contrast(preceeding_trials)==contrast)
            continue
        end
        ops_local.trials = triggers{iF}(iRep)+[-ops_local.n_preceeding:9];
        [~,t1,t2]=fileparts(filenames{iF});
        sn = strcat(t1,t2);
        shift_path = fullfile(shift_data_path,sn);
        data_out{iRep}=compare_speedShift_gainShift(data,ops_local,shift_path);
        %data_out{iRep}=compare_speedShift_gainShift(data,ops_local);
        test_trials = triggers{iF}(iRep)+[4:13];
        ops_local.trials=test_trials;
        [correlation_shifted,correlation_noshift]=test_alignement(data,ops_local,data_out{iRep}.factors,data_out{iRep}.CID);
        data_out{iRep}.correlation_shifted = correlation_shifted;
        data_out{iRep}.correlation_noshift = correlation_noshift;
        
    end
    output{iF}=data_out;
end
%% across sites
allM=[];
allS=[];
allMHat = [];
allSHat = [];
region = 'VISp';
DEPTH =[] ;
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    
    for iRep =1:2
        if isempty(output{iF}{iRep})
            continue
        end
        
        data_out = output{iF}{iRep};
        dcorr = data_out.correlation_shifted - data_out.correlation_noshift;
        idx_region  = startsWith(data_out.region,region);
        idx = data_out.similarity>.4 & startsWith(data_out.region,region)';% & dcorr'>0;
        if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.2

            tmp = squeeze(nanmean(data_out.corrMat(idx,:,:)));  
            tmpS = squeeze(nanmean(data_out.shiftMat(idx,:,:)));
            allM=cat(3,allM,tmp);
            allS = cat(3,allS,tmpS);
            allMHat = cat(3,allMHat,squeeze(nanmean(data_out.corrMatHat(idx,:,:))));
            allSHat = cat(3,allSHat,squeeze(nanmean(data_out.shiftMatHat(idx,:,:))));
            
        end
    end
end
figure
trials2show = 18+[-5:10];
subplot(2,1,1)
imagesc(nanmean(allS(trials2show,trials2show,:),3),[-4 4])
axis image
subplot(2,1,2)
imagesc(nanmean(allSHat(trials2show,trials2show,:),3),[-4 4])
axis image

%%
allShifts = allS(trials2show,trials2show,:);
allShiftsCorrected = allSHat(trials2show,trials2show,:);
shifts_per_cell = mean(nanmean(allShifts(1:6,7:10,:),1),2);
shifts_per_cell_corrected = mean(nanmean(allShiftsCorrected(1:6,7:10,:),1),2);

figure
scatter(shifts_per_cell,shifts_per_cell_corrected,15,'filled')
axis image
grid on
hold on
plot([-15 15],[-15 15],'k')
xlim([-15 15])
ylim([-15 15])
xlabel('shift')
ylabel('speed adjusted shift') 
%% across cells
allM=[];
allS=[];
allMHat = [];
allSHat = [];
DEPTH = [];
FACT = [];
region = 'VISp';
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    
    for iRep =1:2
        if isempty(output{iF}{iRep})
            continue
        end
        
        data_out = output{iF}{iRep};
        dcorr = data_out.correlation_shifted - data_out.correlation_noshift;
        idx_region  = startsWith(data_out.region,region);
        idx = data_out.similarity>0.4 & startsWith(data_out.region,region)'& dcorr'>0.0;% & data_out.factors'<0;
        if nnz(idx)<2
            continue
        end

            tmp = squeeze((data_out.corrMat(idx,:,:)));  
            tmpS = squeeze((data_out.shiftMat(idx,:,:)));
            allM=cat(1,allM,tmp);
            allS = cat(1,allS,tmpS);
            allMHat = cat(1,allMHat,squeeze((data_out.corrMatHat(idx,:,:))));
            allSHat = cat(1,allSHat,squeeze((data_out.shiftMatHat(idx,:,:))));
            DEPTH = cat(1,DEPTH,data_out.depth(idx)');
            FACT = cat(1,FACT,data_out.factors(idx)');
        
    end
end

figure
subplot(2,2,1)
imagesc(squeeze(nanmean(allM(:,trials2show,trials2show),1)),[.3 .75])
axis image
subplot(2,2,2)
imagesc(squeeze(nanmean(allMHat(:,trials2show,trials2show),1)),[.3 .75])
axis image

subplot(2,2,3)
imagesc(squeeze(nanmean(allS(:,trials2show,trials2show),1)),[-3 3])
axis image
subplot(2,2,4)
imagesc(squeeze(nanmean(allSHat(:,trials2show,trials2show),1)),[-3 3])
axis image
%%
allShifts = allS(:,trials2show,trials2show);
allShiftsCorrected = allSHat(:,trials2show,trials2show);
shifts_per_cell = mean(nanmean(allShifts(:,1:6,7:10),3),2);
shifts_per_cell_corrected = mean(nanmean(allShiftsCorrected(:,1:6,7:10),3),2);

cl=[-.55 0.1];
figure
subplot(1,2,1)
scatter(shifts_per_cell,DEPTH,15,FACT,'filled')
set(gca,'YDir','reverse')
ylabel('Dist from surface')
xlabel('shift relative to baseline')
colorbar
xl = get(gca,'XLim');
xlim(round(xl));
set(gca,'CLim',cl)
grid on
subplot(1,2,2)
scatter(shifts_per_cell_corrected,DEPTH,15,FACT,'filled')
set(gca,'YDir','reverse')
ylabel('Dist from surface')
xlabel('shift relative to baseline')
colorbar

xlim(round(xl))
grid on
set(gca,'CLim',cl)
%%
figure
scatter(shifts_per_cell,shifts_per_cell_corrected,15,FACT,'filled')
axis image
grid on
hold on
plot([-15 15],[-15 15],'k')
xlim([-15 15])
ylim([-15 15])
xlabel('shift')
ylabel('speed adjusted shift')  
set(gca,'CLim',[-.3 .3])
colorbar
%%
figure
subplot(1,2,1)
plot([shifts_per_cell,shifts_per_cell_corrected]',[DEPTH,DEPTH]','k')
hold on
scatter(shifts_per_cell,DEPTH,15,FACT,'filled')
set(gca,'YDir','reverse')
xlim(xl)
grid on
set(gca,'CLim',cl)


subplot(1,2,2)
plot([shifts_per_cell,shifts_per_cell_corrected]',[DEPTH,DEPTH]','k')
hold on
scatter(shifts_per_cell_corrected,DEPTH,15,FACT,'filled')
set(gca,'YDir','reverse')
xlim(xl)
grid on
set(gca,'CLim',cl)


