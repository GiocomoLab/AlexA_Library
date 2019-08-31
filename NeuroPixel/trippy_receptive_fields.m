[~,~,time_id]=histcounts(sp.st,post);
stim_id = stim_id+1;
stim_id = round(stim_id);
figure
plot(stim_id)
mov_onsets = strfind(stim_id'>0.5,[0 0 1 1])+2;
mov_offsets = strfind(stim_id'>0.5,[1 1 0 0])+2;
if numel(mov_onsets) > numel(mov_offsets)
    mov_onsets=mov_onsets(1:numel(mov_offsets));
end
hold on
plot(mov_onsets,stim_id(mov_onsets),'ro')
%%
frame_number = zeros(size(stim_id));
for iM = 1:numel(mov_offsets)
    idx = mov_onsets(iM):mov_offsets(iM);
    time_vec = post(idx)-post(idx(1));
    frame_number(idx)=round(time_vec*60+1);
end
    

stim_per_spike = stim_id(time_id);
frame_per_spike = frame_number(time_id);
%% load all movs
movMat = [];
root = 'Y:\giocomo\attialex\stim_movies';
for iM = 0:39
    fn = sprintf('Trippy_test_%d.npy',iM);
    dat = readNPY(fullfile(root,fn));
    movMat=cat(1,movMat,dat);
end
% 
% frs = randsample(1:900,200);
% movs = randsample(1:11,200,true);
% 
% ind = frs + (movs-1)*900;
%%
mm=mean(movMat,1);

good_cells = sp.cids(sp.cgs==2);
for iC = 1:15
    fig_i=0;
    figure
for offset = 0:-1:-7
%offset = -2;
    fig_i=fig_i+1;
    idx = sp.clu==good_cells(iC);
    frames_this = frame_per_spike(idx)+offset;
    mov_this = stim_per_spike(idx);
    mov_this = mov_this -1;
    mov_this(mov_this<0)=0;
    valid_idx = frames_this>0;


        ind = frames_this+(mov_this)*900;
        ind = ind(valid_idx);

        subplot(8,1,fig_i);
        mm_this = mean(movMat(ind,:,:),1)-mm;
        imagesc(squeeze(mm_this),[-5 5])
    %     figure
    %     histogram(ind,1000)
    end
    end
    %%
    pp=shiftdim(movMat,1);
