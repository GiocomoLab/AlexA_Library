chunksize=199; %in bins,so thats 200 cm
stride_start = 1;
binsize=2;
stride = 10;
startVec = stride_start:stride:(200-chunksize+1);
chunksPerTrials = numel(startVec);
region = 'MEC';
gain_to_look_at = .6;
contrast = 100;
[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/users/attialex/Desktop/data');
%%
%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'F:\NP_DATA');
%
p=gcp('nocreate');
if isempty(p)
    parpool(12);
end
%%
savepath_root = '/oak/stanford/groups/giocomo/attialex/Images/xcorrv_range_gains';
%savepath_root = '/users/attialex/tmp/';
savepath = fullfile(savepath_root,sprintf('%s_%.2f_%d',region,gain_to_look_at,contrast));
if ~isfolder(savepath)
    mkdir(savepath)
end

%%

nFiles = numel(filenames);
n_chunks = 1;
template_trials = 1:6;
PEAKS=nan(212,1,nFiles);
SHIFTS = PEAKS;
GAINS = nan(212,nFiles);
MouseID = cell(nFiles,1);
NUnits = nan(nFiles,2);
XCORRS = nan(5,21,nFiles);

%cntr = 0;
parfor iF = 1:nFiles
    
    data = load(filenames{iF});
    template_trials = find(data.trial_gain == .5)';
    trials = 1:max(data.trial);
    nT = numel(trials);
    try
        [peak,shift,~,n_units,~,~,~,xcorrs]=calculatePeakShiftSession(data,trials,chunksize,stride_start,stride,region,0.2,binsize,template_trials);
        
        [~,session_name,~] = fileparts(filenames{iF});
        
        TMP=nan(5,21);
        gains = [0.5 0.6 0.7 0.8 1];
        
        for iG=1:numel(gains)
            tidx = data.trial_gain == gains(iG);
            tmp = nanmean(xcorrs(tidx,:));
            TMP(iG,:)=tmp;
        end
        XCORRS(:,:,iF)=TMP;
        gains = data.trial_gain;
        if nT>212
            peak=peak(1:212,:);
            shift=shift(1:212,:);
            gains = gains(1:212,:);
        end
        if nT<212
            tmpP = nan(212,size(peak,2));
            tmpS = nan(212,size(shift,2));
            tmpG=nan(212,1);
            tmpP(1:nT,:)=peak;
            tmpS(1:nT,:)=shift;
            tmpG(1:nT)=gains;
            peak = tmpP;
            shift=tmpS;
            gains = tmpG;
        end
        PEAKS(:,:,iF)=peak;
        SHIFTS(:,:,iF)=shift;
        GAINS(:,iF)=gains;
        MouseID{iF}=session_name;
        NUnits(iF,:)=n_units;
        
        x=startVec+chunksize/2;
        x = x-1;
        x = x*2;
        fig = figure('visible','off');
        hold on
        gains = [1 0.8 0.7 0.6 0.5];
        cmap = [0 0 0;cool(4)];
        for iG = 1:numel(gains)
            if numel(data.trial_gain)>212
            tidx = find(data.trial_gain(1:212) == gains(iG));
            else
                tidx = find(data.trial_gain == gains(iG));
            end
            %%%
            
            %needs fixing here with the index
            %%%
            plot(tidx,peak(tidx,:),'.','Color',cmap(iG,:))
            
        end
        ylim([0., 0.8])
        xlim([0 numel(peak)])
        
        
        
        drawnow
        saveas(fig,fullfile(savepath,sprintf('%s_%s_%.1f_%d_%d.png',session_name,region,gain_to_look_at,contrast,iF)))
        
        close(fig)
    catch ME
        sprintf('%s: %d',filenames{iF},iF)
        rethrow(ME)
    end
    %end
end
%%


output=struct();
output.GAINS=GAINS;
output.SHIFTS = SHIFTS;
output.PEAKS = PEAKS;
output.XCORRS = XCORRS;
output.region = region;
output.NUnits = NUnits;
%figure
%plot(output.X-4000,nanmean(output.Y))
%hold on
save(fullfile(savepath_root,sprintf('allData_%s_%.1f_%d.mat',region,gain_to_look_at,contrast)),'output','-v7.3')
%%

