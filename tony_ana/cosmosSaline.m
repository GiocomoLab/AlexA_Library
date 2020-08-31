%%
data_folder = '/Users/attialex/Desktop/alex';
cell_regions=load(fullfile(data_folder,'20190521-cux2m4293_saline_day2_1'),'cell_regions');
cell_regions = cell_regions.cell_regions;
trial_types=load(fullfile(data_folder,'20190521-cux2m4293_saline_day2_1'),'trial_types');
trial_types = trial_types.trial_types;
cell_hemis=load(fullfile(data_folder,'20190521-cux2m4293_saline_day2_1'),'cell_hemis');
cell_hemis =cell_hemis.cell_hemis;
act = load(fullfile(data_folder,'20190521-cux2m4293_saline_day2_1'),'trial_C');
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
respPreSaline={};
respPostSaline={};
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
    respPreSaline{jj} = mean(snps(:,:,1:6),3); % for each neuron, average across all trials
    respPostSaline{jj} = mean(snps(:,:,13:18),3); % for each neuron, average across all trials
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
for iStim = 1:length(respPreSaline)
    
    for iR = 5
        IDX = ismember(cell_regions,regions{iR});
        subplot(3,3,iStim)
                params = struct();
        params.winIDX = (1:size(respPreSaline{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate; % convert sample number to time
    params.xLim = [-1 4]; % shows only from -1s to 3s
        plotAVGSEM(respPreSaline{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89)
        hold on
        params.winIDX = (1:size(respPreSaline{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate;
        plotAVGSEM(respPostSaline{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89,'col',[1 0 0])
        title(stimNames{iStim})
    end
end
%% vis, average gratings
nFramesPre = Inf;
nFramesPost = Inf;
for iStim = 5:8
    nFramesPre=min(nFramesPre,size(respPreSaline{iStim},2));
    nFramesPost = min(nFramesPost,size(respPreSaline{iStim},2));
end
GratSalinePre = [];
GratSalinePost=[];
for iStim = 5:8
    GratSalinePre = cat(3,GratSalinePre,respPreSaline{iStim}(:,1:nFramesPre));
    GratSalinePost = cat(3,GratSalinePost,respPostSaline{iStim}(:,1:nFramesPost));
end
%%
GratAmpPreSaline = squeeze(mean(GratSalinePre(:,95:155,:),2)-mean(GratSalinePre(:,39:89,:),2));
GratAmpPostSaline = squeeze(mean(GratSalinePost(:,95:155,:),2)-mean(GratSalinePost(:,39:89),2));
thresh = 10;
responsiveIDX = any(GratAmpPreSaline>thresh,2) & any(GratAmpPostSaline>thresh,2);
%%
for iGrat = 1:4
    IDX = ismember(cell_regions,'VIS');
    tmp = GratSalinePre(:,:,iGrat);
    tmp = bsxfun(@minus,tmp,GratAmpPreSaline(:,iGrat));
    tmp = tmp(IDX,:);
    figure
    [~,sidx]=sort(GratAmpPreSaline(IDX,iGrat),'descend');
    imagesc(tmp(sidx,:),[-100 100]);
end
%% vis, all stim, positive grating response pre and post
regions = unique(cell_regions);
figure
for iStim = 1:length(respPreSaline)
    
    for iR = 5
        IDX = ismember(cell_regions,regions{iR}) & responsiveIDX';
        subplot(3,3,iStim)
                params = struct();
        params.winIDX = (1:size(respPreSaline{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate; % convert sample number to time
    params.xLim = [-1 4]; % shows only from -1s to 3s
        plotAVGSEM(respPreSaline{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89)
        hold on
        params.winIDX = (1:size(respPreSaline{iStim},2))-89; % in sample number, has to be the same as number of frames in resp
    params.masterTime = params.winIDX/frameRate;
        plotAVGSEM(respPostSaline{iStim}(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89,'col',[1 0 0])
        title(stimNames{iStim})
    end
end
%%
salineData.GratAmpPreSaline=GratAmpPreSaline;
salineData.GratAmpPostSaline=GratAmpPostSaline;
salineData.cell_hemis = cell_hemis;
salineData.cell_regions=cell_regions;
salineData.GratSalinePre = GratSalinePre;
salineData.GratSalinePost = GratSalinePost;
