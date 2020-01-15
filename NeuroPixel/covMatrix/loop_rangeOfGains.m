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
        catch ME
            sprintf('%s: %d',filenames{iF},iF)
            rethrow(ME)
        end
        [~,session_name,~] = fileparts(filenames{iF});
        
        gains = [0.5 0.6 0.7 0.8 1];
        for iG=1:numel(gains)
            tidx = data.trial_gain == iG;
            tmp = nanmean(xcorrs(tidx,:));
            XCORRS(iG,:,iF)=tmp;
        end
        
        if nT>212
            peak=peak(1:212,:);
            shift=shift(1:212,:);
        end
        if nT<212
            tmpP = nan(nT,size(peak,2));
            tmpS = nan(nT,size(shift,2));
            tmpP(1:nT,:)=peak;
            tmpS(1:nT,:)=shift;
            peak = tmpP;
            shift=tmpS;
        end
        PEAKS(:,:,iF)=peak;
        SHIFTS(:,:,iF)=shift;
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

                tidx = find(data.trial_gain == gains(iG));
                plot(tidx,peak(tidx,:),'.','Color',cmap(iG,:))
            
        end
        ylim([0., 0.8])
        xlim([0 numel(peak)])
      
        
        
        drawnow
        saveas(fig,fullfile(savepath,sprintf('%s_%s_%.1f_%d_%d.png',session_name,region,gain_to_look_at,contrast,iF)))
        
        close(fig)
    %end
end
%%
X=[];


x=startVec+chunksize/2;
x = x-1;
x = x*binsize;
for iT = 1:numel(tt)
    tmp = x+400*(iT-1);
    X=cat(2,X,tmp);
    %plot(tmp,peak(iT,:),'b.')
end
Y=zeros(size(PEAKS,3),numel(X));
S = Y;
chunks_per_trial = size(PEAKS,2);
for iS = 1:size(PEAKS,3)
    for iT = 1:numel(tt)
        idx = ((iT-1)*chunks_per_trial+1):iT*chunks_per_trial;
        tmp = squeeze(PEAKS(iT,:,iS));
        Y(iS,idx)=tmp;
        tmp = squeeze(SHIFTS(iT,:,iS));
        
        S(iS,idx)=tmp;
    end
end

output=struct();
output.X=X;
output.Y = Y;
output.S = S;
output.region = region;
output.gain = gain_to_look_at;
output.contrast = contrast;
output.loop_data = loop_data;
output.NUnits = NUnits;
output.XTX =nanmean(XTX,3);
output.YYT = YYT;
output.SPEED = SPEED;
output.FR = FR;
%figure
%plot(output.X-4000,nanmean(output.Y))
%hold on
save(fullfile(savepath_root,sprintf('allData_%s_%.1f_%d.mat',region,gain_to_look_at,contrast)),'output','-v7.3')
%%

