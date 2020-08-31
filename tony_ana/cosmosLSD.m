%%
data_folder = '/Users/attialex/Desktop/alex';
filename = '20190522-cux2m4293_lsd_1-003';
cell_regions=load(fullfile(data_folder,filename),'cell_regions');
cell_regions = cell_regions.cell_regions;
trial_types=load(fullfile(data_folder,filename),'trial_types');
trial_types = trial_types.trial_types;
cell_hemis=load(fullfile(data_folder,filename),'cell_hemis');
cell_hemis =cell_hemis.cell_hemis;
act = load(fullfile(data_folder,filename),'trial_C');
fn=fieldnames(act);
act=act.(fn{1});
%%
nCells = length(cell_regions);
RSP_cells = ismember(cell_regions,'RSP');
MO_cells = ismember(cell_regions,'MO');
PTLp_cells = ismember(cell_regions,'PTLp');
SSp_cells = ismember(cell_regions,'SSp');
stimNamesTemp = unique(trial_types);
stimNames = {stimNamesTemp{1:8},stimNamesTemp{10}};
regionNames = unique(cell_regions);

frameRate = 1/0.034;
col = [0, 0, 1; 1 0 1; 1 0 0];
figure('color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
respPrelsd={};
respPostlsd={};
for jj = 1:length(stimNames)
    
    stimIdx =  find(ismember(trial_types,stimNames{jj}));
    region_cells = true(size(cell_hemis'));%ismember(cell_regions,'VIS') & logical(cell_hemis');
    % find minimal window length for specific stimulus
    window_length_temp = [];
    for trialIdx = 1:length(stimIdx)
        window_length_temp(trialIdx) = size(act{stimIdx(trialIdx)},2); %#ok<SAGROW>
    end
    window_length = min(window_length_temp);
    
    % get snps data (neurons x window length x trials)
    snps = zeros(nCells,window_length+round(3*frameRate),length(stimIdx));
    for trialIdx = 1:length(stimIdx)
        stim__pregrey = act{stimIdx(trialIdx)-1}(:,size(act{stimIdx(trialIdx)-1},2)-round(3*frameRate)+1:end);
        stim__postgrey = act{stimIdx(trialIdx)}(:,1:window_length);
        snps(:,:,trialIdx) = cat(2,stim__pregrey, stim__postgrey);
    end
    respPrelsd{jj} = mean(snps(:,:,1:6),3); % for each neuron, average across all trials
    respPostlsd{jj} = mean(snps(:,:,13:18),3); % for each neuron, average across all trials
end
%% all regions, all stim
% regions = unique(cell_regions);
% for iStim = 1:length(respPre)
%     figure
%     for iR = 1:numel(regions)
%         IDX = ismember(cell_regions,regions{iR});
%         subplot(1,numel(regions),iR)
%         plot(mean(respPre{iStim}(IDX,:)))
%         hold on
%         plot(mean(respPost{iStim}(IDX,:)))
%         title(regions{iR})
%     end
% end

%% vis, all stim
regions = unique(cell_regions);
figure
for iStim = 1:length(respPrelsd)
    
    for iR = 5
        IDX = ismember(cell_regions,regions{iR});
        subplot(3,3,iStim)
                params = struct();
        params.winIDX = (1:size(respPrelsd{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate; % convert sample number to time
    params.xLim = [-1 4]; % shows only from -1s to 3s
        plotAVGSEM(respPrelsd{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89)
        hold on
        params.winIDX = (1:size(respPrelsd{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate;
        plotAVGSEM(respPostlsd{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89,'col',[1 0 0])
        title(stimNames{iStim})
    end
end
%% vis, average gratings
nFramesPre = Inf;
nFramesPost = Inf;
for iStim = 5:8
    nFramesPre=min(nFramesPre,size(respPrelsd{iStim},2));
    nFramesPost = min(nFramesPost,size(respPrelsd{iStim},2));
end
GratlsdPre = [];
GratlsdPost=[];
for iStim = 5:8
    GratlsdPre = cat(3,GratlsdPre,respPrelsd{iStim}(:,1:nFramesPre));
    GratlsdPost = cat(3,GratlsdPost,respPostlsd{iStim}(:,1:nFramesPost));
end
%%
GratAmpPrelsd = squeeze(mean(GratlsdPre(:,95:155,:),2)-mean(GratlsdPre(:,39:89,:),2));
GratAmpPostlsd = squeeze(mean(GratlsdPost(:,95:155,:),2)-mean(GratlsdPost(:,39:89),2));
thresh = 10;
responsiveIDX = any(GratAmpPrelsd>thresh,2) & any(GratAmpPostlsd>thresh,2);
%%
for iGrat = 1:4
    IDX = ismember(cell_regions,'VIS');
    tmp = GratlsdPre(:,:,iGrat);
    tmp = bsxfun(@minus,tmp,GratAmpPrelsd(:,iGrat));
    tmp = tmp(IDX,:);
    figure
    [~,sidx]=sort(GratAmpPrelsd(IDX,iGrat),'descend');
    imagesc(tmp(sidx,:),[-100 100]);
end
%% vis, all stim, positive grating response pre and post
regions = unique(cell_regions);
figure
for iStim = 1:length(respPrelsd)
    
    for iR = 5
        IDX = ismember(cell_regions,regions{iR}) & responsiveIDX';
        subplot(3,3,iStim)
                params = struct();
        params.winIDX = (1:size(respPrelsd{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate; % convert sample number to time
    params.xLim = [-1 4]; % shows only from -1s to 3s
        plotAVGSEM(respPrelsd{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89)
        hold on
        params.winIDX = (1:size(respPrelsd{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate;
        plotAVGSEM(respPostlsd{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89,'col',[1 0 0])
        title(stimNames{iStim})
    end
end
%%
lsdData.GratAmpPrelsd=GratAmpPrelsd;
lsdData.GratAmpPostlsd=GratAmpPostlsd;
lsdData.cell_hemis = cell_hemis;
lsdData.cell_regions=cell_regions;
lsdData.GratlsdPre = GratlsdPre;
lsdData.GratlsdPost = GratlsdPost;

