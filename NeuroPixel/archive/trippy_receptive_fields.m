[~,~,time_id]=histcounts(data_trippy.sp.st,data_trippy.vr_data.time); %for each spike, which time bin it is in
stim_id = data_trippy.vr_data.stimID+1; % which movie is displayed
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
post=data_trippy.vr_data.time;
frame_number = zeros(size(stim_id));
for iM = 1:numel(mov_offsets)
    idx = mov_onsets(iM):mov_offsets(iM);
    time_vec = post(idx)-post(idx(1));
    frame_number(idx)=round(time_vec*60+1);
end


stim_per_spike = stim_id(time_id); %for each spike, which stimulus movie was displayed
frame_per_spike = frame_number(time_id); % for each spike, which frame number was displayed
%% load all movs
movMat = [];
movCell = {};
root = 'F:\Alex\stim_movies';
for iM = 0:39
    fn = sprintf('Trippy_test_%d.npy',iM);
    dat = readNPY(fullfile(root,fn));
    movMat=cat(1,movMat,dat);
    movCell{iM+1}=shiftdim(dat,1);
end
%
% frs = randsample(1:900,200);
% movs = randsample(1:11,200,true);
%
% ind = frs + (movs-1)*900;
%%
mov2 = ones(96,176,length(post),'uint8')*128;
% t_mov =
% for iR = 1:96
%     for iC=1:176
%         v = interp1(t_mov,movMat(iR,iC),xq);
%         mov2(:,iR,iC)=v;
%     end
% end
%
onset_number = zeros(length(post),1);
for iT=1:numel(stim_id)
    mov_id = stim_id(iT);
    id = find(iT>=mov_onsets & iT<=mov_offsets);
    if ~isempty(id) & mov_id
        onset_number(iT)=id;
        
        start_this = mov_onsets(id);
        stop_this = mov_offsets(id);
        dd = (stop_this-start_this);
        elapsed = iT-start_this;
        if elapsed==0
            elapsed = 1;
        end
        onset_number(iT) = ceil(elapsed/dd*900);
        mov2(:,:,iT)=movCell{mov_id}(:,:,onset_number(iT));
    end
    %frame_number =round((stop_this-start_this)/900*elapsed);
end
%%
good_cells = data_spatial.sp.cids(data_spatial.sp.cgs==2 & startsWith(data_spatial.anatomy.cluster_parent,'VISp'));

[fr,spikeCount] = calcFRVsTime(good_cells,data_trippy,load_mismatch_opt);
%%
spikeCountT=spikeCount';
CM = {};
k = [0.2:0.02:0.4];
shift = round(k*50);
for iS=1:numel(shift)
    corrMat = zeros(size(mov2,1),size(mov2,2),numel(good_cells));
    k=shift(iS);
    mov2D=double(mov2);
    for iR=1:2:size(mov2,1)
        for iC=1:2:size(mov2,2)
            pp=squeeze(mov2D(iR,iC,:));
            pp = circshift(pp,k,1);
            corrMat(iR,iC,:)=corr(double(pp),spikeCountT);
        end
    end
    CM{iS}=corrMat;
end
%%
figure
corrMat = CM{9};
for ii=1:104
    imagesc(squeeze(corrMat(1:2:96,1:2:176,ii)),[-.01 .1])
    pause
    cla
end
%%
figure
for iC=1:6%numel(good_cells)
    figure
    for ii=1:11
        corrMat = CM{ii};
        subplot(3,4,ii)
        imagesc(squeeze(corrMat(1:2:96,1:2:176,iC)),[-.01 .01])
    end
    %pause
    %clf
end
%%
opt = load_mismatch_opt;
opt.max_lag = 30;
[cm,frMat]=trialCorrMat(good_cells,5:20,data_spatial,opt);
tc=squeeze(mean(frMat,2));

%%
visp_clu = data_spatial.sp.cids(data_spatial.sp.cgs==2 & startsWith(data_spatial.anatomy.cluster_parent,'VISp'));
%good_cells = sp.cids(sp.cgs==2 & startsWith(data.anatomy.parent_shifted,'VISp'));
good_cells = visp_clu;
mm=squeeze(mean(mov2,3));
seed_neuron = 6;
pwd = squareform(pdist(tc,'Cosine'));
dist_to_seed = pwd(seed_neuron,:);
[a,sid]=sort(dist_to_seed);
offsets = 0.1:0.016:0.3;
corrMat = zeros(size(mov2,1),size(mov2,2),numel(offsets),numel(good_cells));

for iC = 1:5
    fig_i=0;
    figure
    c_idx = sid(iC);
    
    for offset = offsets
        %offset = -2;
        fig_i=fig_i+1;
        idx = data_trippy.sp.clu==good_cells(c_idx);
        st = data_trippy.sp.st(idx)-offset;
        v_idx =st>=0;
        [~,~,time_idx]=histcounts(st(v_idx),post);
        
        
        subplot(3,5,fig_i);
        mm_this = mean(mov2(:,:,time_idx),3)-mm;
        imagesc(squeeze(mm_this),[-5 5])
        
    end
    subplot(3,5,15)
    plot(tc(c_idx,:))
end
%% all RFS with spike triggred average
mm=squeeze(mean(mov2,3));
offsets = 0.1:0.016:0.3;
staMat = zeros(size(mov2,1),size(mov2,2),numel(offsets),numel(good_cells));

for iC = 1:numel(good_cells)
    
    iO=0;
    for offset = offsets
        iO=iO+1;
        idx = data_trippy.sp.clu==good_cells(iC);
        st = data_trippy.sp.st(idx)-offset;
        v_idx =st>=0;
        [~,~,time_idx]=histcounts(st(v_idx),post);
        
        
        mm_this = mean(mov2(:,:,time_idx),3)-mm;
        staMat(:,:,iO,iC)=mm_this;
    end
    
end
%%
rf = mean(staMat(:,:,2:3,:),3);
figure
for ii=1:numel(good_cells)
    for iD=1:size(staMat,3)
        subplot(3,5,iD)
        imagesc(squeeze(staMat(:,:,iD,ii)),[-5 5])
    end
    subplot(3,5,14)
    %tmp = squeeze(rf(:,:,:,ii));
    %tmp = 0.5*imgaussfilt(squeeze(staMat(:,:,2,ii)))+0.5*imgaussfilt(squeeze(staMat(:,:,3,ii)));
    %ff=imgaussfilt3(squeeze(staMat(:,:,3:7,ii)),[0.5,0.5,00.1]);
    ff= squeeze(staMat(:,:,3:7,ii));
    tmp = mean(ff,3);
    tmp = imgaussfilt(tmp);
    Z = zscore(tmp,[],'All');
    Z(abs(Z)<2.5)=0;
    
    imagesc(Z,[-5 5])
    hold on
    p = fit_receptive_field(tmp,2.5);
    for iR=1:numel(p)
        if p{iR}.field_sign ==-1
            col='b';
        else
            col = 'r';
        end
        plot(p{iR}.xy(:,1),p{iR}.xy(:,2),col)
    end
    pause
    cla
end
%%
rf =zeros(size(staMat,1),size(staMat,2),size(staMat,3));
figure
hold on
for ii=1:numel(good_cells)
   
    ff= squeeze(staMat(:,:,3:7,ii));
    tmp = mean(ff,3);
    tmp = imgaussfilt(tmp);
    Z = zscore(tmp,[],'All');
    Z(abs(Z)<2.5)=0;
    rf(:,:,ii)=Z;
    
     p = fit_receptive_field(tmp,2.5);
    for iR=1:numel(p)
        if p{iR}.field_sign ==-1
            col='b';
        else
            col = 'r';
        end
        subplot(1,2,1)
        hold on
        plot(p{iR}.xy(:,1),p{iR}.xy(:,2),col)
        subplot(1,2,2)
        hold on
        plot(p{iR}.mu(1),p{iR}.mu(2),strcat('x',col'))
    end
end
%%
gg=reshape(rf,[],size(rf,3));
pwd = (pdist(gg','Cosine'));
%%
[~,sid]=sort(pwd);
[ii,jj]=triind2sub([103,103],sid);
%%

for iC=1:10
    figure
    subplot(2,2,1)
    imagesc(rf(:,:,ii(iC)),[-5 5])
    subplot(2,2,2)
    imagesc(rf(:,:,jj(iC)),[-5 5])
    subplot(2,2,3)
    plot(tc(ii(iC),:))
    hold on
    plot(tc(jj(iC),:))
end

%%


mm=mean(movMat,1);

good_cells = sp.cids(sp.cgs==2 & startsWith(data.anatomy.parent_shifted,'VISp'));
for iC = 1:15
    fig_i=0;
    figure
    for offset = -5:-3:-15
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
