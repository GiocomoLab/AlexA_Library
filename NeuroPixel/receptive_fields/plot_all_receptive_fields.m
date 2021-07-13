data_table = readtable("C:\Users\giocomolab\Downloads\vi_trippy_pairs.xlsx");
trippy_path = 'F:\Alex\new_2';
v1_path = 'Z:\giocomo\attialex\NP_DATA_corrected';
im_save_root = 'F:\Alex\receptive_fields_2';
%%
%% load all movs
% movMat = [];
% movCell = {};
% root = 'F:\Alex\stim_movies';
% for iM = 0:39
%     fn = sprintf('Trippy_test_%d.npy',iM);
%     dat = readNPY(fullfile(root,fn));
%     movMat=cat(1,movMat,dat);
%     movCell{iM+1}=shiftdim(dat,1);
% end
mov = load('F:\Alex\stim_movies\stim_movies.mat');
movCell = mov.movCell;
%%
cm = brewermap(20,'*RdBu');
for iR = 1:size(data_table,1)
    try
    data_trippy = load(fullfile(trippy_path,data_table.Trippy_name{iR}));
    
    data_spatial = load(fullfile(v1_path,data_table.Spatial_name{iR}));
    im_save_path = fullfile(im_save_root,data_table.Trippy_name{iR});
    if ~isfolder(im_save_path)
        mkdir(im_save_path);
    end
    
%     if isfield(data_spatial.anatomy,'parent_shifted')
%         reg = data_spatial.anatomy.parent_shifted;
%     else
%         reg = data_spatial.anatomy.cluster_parent;
%     end
%     visp_cells = startsWith(reg,'VISp') & data_spatial.sp.cgs==2;
    
%    good_cells = data_spatial.sp.cids(visp_cells);
    good_cells = data_spatial.sp.cids(data_spatial.sp.cgs==2);
    
    [~,~,time_id]=histcounts(data_trippy.sp.st,data_trippy.vr_data.time); %for each spike, which time bin it is in
    stim_id = data_trippy.vr_data.stimID+1; % which movie is displayed
    stim_id = round(stim_id);
    
    mov_onsets = strfind(stim_id'>0.5,[0 0 1 1])+2;
    mov_offsets = strfind(stim_id'>0.5,[1 1 0 0])+2;
    if numel(mov_onsets) > numel(mov_offsets)
        mov_onsets=mov_onsets(1:numel(mov_offsets));
    end
    
    %calculate movie for this rep
    post=data_trippy.vr_data.time;
    
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
        if ~isempty(id) && mov_id
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
    
    mm=squeeze(mean(mov2,3));
    offsets = 0.1:0.016:0.3;
    staMat = zeros(size(mov2,1),size(mov2,2),numel(offsets),numel(good_cells));
    fig = figure('Position',[ 172         277        1442         701]);
    fields = cell(numel(good_cells),1);
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
        
        for iD=1:size(staMat,3)
            subplot(3,5,iD)
            imagesc(squeeze(staMat(:,:,iD,iC)),[-5 5])
            colormap(cm)
            axis image
            

        end
        subplot(3,5,14)
        %tmp = squeeze(rf(:,:,:,ii));
        %tmp = 0.5*imgaussfilt(squeeze(staMat(:,:,2,ii)))+0.5*imgaussfilt(squeeze(staMat(:,:,3,ii)));
        %ff=imgaussfilt3(squeeze(staMat(:,:,3:7,ii)),[0.5,0.5,00.1]);
        ff= squeeze(staMat(:,:,3:7,iC));
        tmp = mean(ff,3);
        tmp = imgaussfilt(tmp);
        Z = zscore(tmp,[],'All');
        Z(abs(Z)<2.5)=0;
        
        imagesc(Z,[-5 5])
        axis image
        hold on
        p = fit_receptive_field(tmp,2.5);
        for iROI=1:numel(p)
            if p{iROI}.field_sign ==-1
                col='g';
            else
                col = 'k';
            end
            plot(p{iROI}.xy(:,1),p{iROI}.xy(:,2),col)
        end
    fields{iC}=p;
    saveas(fig,fullfile(im_save_path,sprintf('clu_%d.png',good_cells(iC))))
    clf
    end
    close(fig)
    save(fullfile(im_save_path,'receptive_fields.mat'),'staMat','fields','good_cells')
    catch ME
        sprintf('Did not work for %s',data_table.Trippy_name{iR})
    end
end