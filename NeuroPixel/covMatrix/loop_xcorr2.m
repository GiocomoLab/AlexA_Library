chunksize=100; %in bins,so thats 200 cm
stride_start = 10;
binsize=2;
stride = 10;
startVec = stride_start:stride:(200-chunksize+1);
chunksPerTrials = numel(startVec);
region = 'VISp';
gain_to_look_at = [];
contrast = 50
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
savepath_root = '/oak/stanford/groups/giocomo/attialex/Images/xcorrv_CONTRAST2';
%savepath_root = '/users/attialex/tmp/';
savepath = fullfile(savepath_root,sprintf('%s_%.2f_%d',region,gain_to_look_at,contrast));
if ~isfolder(savepath)
    mkdir(savepath)
end

%%
n_chunks = 0;
chunk_idx = triggers;
loop_data = struct();
for ii=1:numel(filenames)
    %n_chunks = n_chunks+numel(triggers{ii});
    if numel(triggers{ii})>4
        triggers{ii} = randsample(triggers{ii},4);
    end
        
    for iC=1:numel(triggers{ii})
        n_chunks = n_chunks+1;
        chunk_idx{ii}(iC)=n_chunks;
        loop_data(n_chunks).filename = filenames{ii};
        loop_data(n_chunks).trigger = triggers{ii}(iC);
    end
end

tt=(-6:9);
template_trials = 1:6;
PEAKS=nan(numel(tt),chunksPerTrials,n_chunks);
SHIFTS = PEAKS;
MouseID = cell(n_chunks,1);
NUnits = nan(n_chunks,2);
XTX = zeros(numel(tt)*200,numel(tt)*200,n_chunks);
YYT = zeros(numel(tt),numel(tt),n_chunks);
SPEED=zeros(numel(tt),200,n_chunks);
%cntr = 0;
parfor iF = 1:n_chunks
%     data = load(filenames{iF});
%     for iTrigger = 1:10
%         cntr = chunk_idx{iF}(iTrigger);
%         current_trig = triggers{iF}(iTrigger);
        data = load(loop_data(iF).filename);
        current_trig = loop_data(iF).trigger;
        
        trials = current_trig+tt;
        peak = 0;
        shift = 0;
        xtx = 0;
        try
        [peak,shift,xtx,n_units,yyt,speed_mat]=calculatePeakShiftSession(data,trials,chunksize,stride_start,stride,region,0.2,binsize,template_trials);
        catch ME
            sprintf('%s: %d',loop_data(iF).filename,iF)
            rethrow(ME)
        end
        [~,session_name,~] = fileparts(loop_data(iF).filename);

        PEAKS(:,:,iF)=peak;
        SHIFTS(:,:,iF)=shift;
        MouseID{iF}=session_name;
        NUnits(iF,:)=n_units;
        XTX(:,:,iF)=xtx;
        YYT(:,:,iF)=yyt;
        SPEED(:,:,iF)=speed_mat;
        x=startVec+chunksize/2;
        x = x-1;
        x = x*2;
        fig = figure('visible','off');
        subplot(1,2,1)
        hold on
        for iT = 1:numel(tt)
            tmp = x+400*(iT-1)-10*400;
            if ismember(iT,[abs(tt(1))+(1:4)])
                plot(tmp,peak(iT,:),'r.')
                
            else
                
                plot(tmp,peak(iT,:),'b.')
            end
        end
        ylim([0., 0.8])
        
        subplot(1,2,2)
        imagesc(xtx,[0 1]);
        xline(abs(min(tt))*400/binsize,'r');
        yline(abs(min(tt))*400/binsize,'r');
        axis image;
        
        
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
%figure
%plot(output.X-4000,nanmean(output.Y))
%hold on
save(fullfile(savepath_root,sprintf('allData_%s_%.1f_%d.mat',region,gain_to_look_at,contrast)),'output','-v7.3')
%%

