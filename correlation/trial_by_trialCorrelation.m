% filenames={'npF2_1015_contrasttrack_gainchanges_2.mat',...
%     'npF2_1016_contrasttrack_gainchanges_1.mat',...
%     'npF3_1018_contrasttrack_gainchanges_1.mat',...
%     'npF3_1019_contrasttrack_gainchanges_contrast_1.mat',...
%     'npF4_1023_gaincontrast_1.mat',...
%     'npF4_1025_gaincontrast_2.mat '};
%restoredefaultpath
if ispc()
    addpath(genpath('C:\code\AlexA_Library'));
    addpath(genpath('C:\code\boundedline'));
    addpath(genpath('F:\code\cortexlab_spikes'));
    
    % filenames = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
    %     'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    %     'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    %     'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    %     'G5/1210_mismatch_1/1210_mismatch_1.mat'
    %     };
    
    filenames = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');
    filenames = dir('Z:\giocomo\attialex\NP_DATA\*mismatch*.mat');
    
    root_dir='F:\';
else
    run('/home/users/attialex/AlexA_Library/default_paths.m')
    filenames=dir(fullfile(OAK,'attialex','NP_DATA','*gain*.mat'));
end
%%
aggregateData = struct();
aggregateData.correlationMatrix={};
aggregateData.trial_gain = {};
aggregateData.trial_contrast = {};
aggregateData.parent = {};

%%
for iF = 1:numel(filenames)
    try
        clear anatomy
        load(fullfile(filenames(iF).folder,filenames(iF).name));
        
        trials=1:max(trial);
        spatialMap=[];
        dwell_time=[];
        edges=[0:5:410];
        edges(1)=-.01;
        posx(posx<0)=0;
        for iT=1:length(trials)
            idxVR=trial==trials(iT);
            t_time=post(idxVR);
            start=min(t_time);
            stop=max(t_time);
            idxNP=sp.st<stop & sp.st>=start;
            [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.clu)+1);
            spatialMap=cat(3,spatialMap,spM);
            dwell_time=cat(1,dwell_time,dT);
        end
        %cellIDX=find(sp.cgs>=1);
        spatialMap=spatialMap(sp.cids+1,:,:);
        spatialMap=spatialMap(sp.cgs==2,:,:);
        spatialMap=spatialMap(:,1:end-1,:);
        dwell_time=dwell_time(:,1:end-1);
        %normalize by dwell time in each bin
        dt=dwell_time';
        dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
        for ii=1:size(spatialMap,1)
            spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
        end
        spatialMap(isnan(spatialMap))=0;
        
        
        
        correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
        diagAll=zeros(size(spatialMap,1),size(spatialMap,3)-1);
        for iC=1:size(spatialMap,1)
            tmp=corr(squeeze(spatialMap(iC,:,:)));
            correlation_All(:,:,iC)=tmp;
            diagAll(iC,:)=diag(tmp,1);
        end
        
        if exist('anatomy','var')
            if isfield(anatomy,'region_shifted')
                tmp_region = anatomy.region_shifted;
                tmp_parent = anatomy.parent_shifted;
            else
                if isfield(anatomy,'cluster_region')
                    tmp_region = anatomy.cluster_region;
                else % because for now we only have cluster parent for mec data
                    tmp_region = anatomy.cluster_parent;
                end
                
                tmp_parent = anatomy.cluster_parent;
            end
            tmp_parent = tmp_parent(sp.cgs==2);
            tmp_region = tmp_region(sp.cgs==2);
            if numel(tmp_region) ~= nnz(sp.cgs==2)
                error('anatomy and real clusters do not match')
            end
        else
            tmp_region = cell(1,nnz(sp.cgs==2));
            tmp_parent = cell(1,nnz(sp.cgs==2));
        end
        if ~isrow(tmp_parent)
            tmp_parent = tmp_parent';
        end
        if ~isrow(tmp_region)
            tmp_region = tmp_region';
        end
        
        aggregateData.correlationMatrix{iF}=correlation_All;
        aggregateData.trial_gain{iF}=trial_gain;
        aggregateData.trial_contrast{iF}=trial_contrast;
        aggregateData.session{iF}=filenames(iF).name;
        aggregateData.parent{iF}=tmp_parent;
    catch ME
        disp(ME.message)
        disp(sprintf('error with file %d, %s',iF,filenames(iF).name))
    end
end
save(fullfile(OAK,'attialex',strcat('TBT_aggregate',date,'.mat')),'aggregateData','-v7.3')

%% all cells
mat=[];
trials = [];
for ii=1:numel(aggregateData.session)
    if ~strfind(aggregateData.session{ii},'gain')
        continue
    end
    if ismember(aggregateData.session{ii},'contrast')
        continue
    end
    if numel(unique(aggregateData.trial_gain{ii}))<5
        continue
    end
    if numel(aggregateData.trial_gain{ii})<212
        continue
    end
    mat=cat(3,mat,nanmean(aggregateData.correlationMatrix{ii}(1:212,1:212,:),3));
end
figure
imagesc(mean(mat,3),[0 .2])
colorbar
axis image

%% all cells
mat=[];
trials = [];
for ii=1:numel(aggregateData.session)
    if ~strfind(aggregateData.session{ii},'gain')
        continue
    end
    if ismember(aggregateData.session{ii},'contrast')
        continue
    end
    if numel(unique(aggregateData.trial_gain{ii}))<5
        continue
    end
    if numel(aggregateData.trial_gain{ii})<212
        continue
    end
    data = aggregateData.correlationMatrix{ii}(1:212,1:212,:);
    stability=[];
    bl_trials = aggregateData.trial_gain{ii}(1:212)==1;
    for iC=1:size(data,3)
        tmp = diag(data(bl_trials,bl_trials,iC),1);
        stability(iC)=mean(tmp);
    end
    valididx = stability>.1;
    siteidx = strcmp(aggregateData.parent{ii},'MEC');
    valididx = siteidx & valididx;
    if nnz(valididx)<3
        continue
    end
    frac = nnz(valididx)/numel(valididx);
    disp(frac);
    mat=cat(3,mat,nanmean(data(:,:,valididx),3));
end
figure
imagesc(mean(mat,3),[0 .4])
colorbar
axis image


