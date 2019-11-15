chunksize=100; %in bins,so thats 200 cm
stride_start = 10;
binsize=2;
stride = 10;
startVec = stride_start:stride:(200-chunksize+1);
chunksPerTrials = numel(startVec);
region = 'MEC';
contrast = 100;
gain_to_look_at = .5;

[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/users/attialex/Desktop/data');

%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'F:\NP_DATA');
%
p=gcp('nocreate');
if isempty(p)
    parpool(8);
end

savepath_root = '/oak/stanford/groups/giocomo/attialex/Images/xcorrv6';
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

tt=(-8:13);

PEAKS=nan(numel(tt),chunksPerTrials,n_chunks);
SHIFTS = PEAKS;
%cntr = 0;
parfor iF = 1:n_chunks
%     data = load(filenames{iF});
%     for iTrigger = 1:10
%         cntr = chunk_idx{iF}(iTrigger);
%         current_trig = triggers{iF}(iTrigger);
        data = load(loop_data(iF).filename);
        current_trig = loop_data(iF).trigger;
        
        trials = current_trig+tt;
        try
        [peak,shift,xtx]=calculatePeakShiftSession(data,trials,chunksize,stride_start,stride,region,0.4,binsize);
        catch ME
            sprintf('%s: %d',loop_data(iF).filename,iF)
            rethrow(ME)
        end
        PEAKS(:,:,iF)=peak;
        SHIFTS(:,:,iF)=shift;
        
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
        
        [~,session_name,~] = fileparts(loop_data(iF).filename);
        
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
figure
plot(output.X-4000,nanmean(output.Y))
hold on
save(fullfile(savepath_root,sprintf('allData_%s_%.1f_%d.mat',region,gain_to_look_at,contrast)),'output')
%%

