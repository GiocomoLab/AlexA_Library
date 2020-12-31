function [glmData] = fitGLM_OLCL(dataOL,dataCL,good_cells,downsample)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% for OL, used to fit model

opt = load_mismatch_opt;
velM = preprocessSignal(dataOL.vr_data_resampled.velM',21);
velP = preprocessSignal(dataOL.vr_data_resampled.velP',21);
pa = removeOutliers(dataOL.pupil_area,[1 99]);
pupil_area = preprocessSignal(pa,11);

p_x = removeOutliers(dataOL.pupil_com(:,1),[1 99]);
pupil_x = preprocessSignal(p_x,11);
p_y = removeOutliers(dataOL.pupil_com(:,2),[1 99]);
pupil_y = preprocessSignal(p_y,11);
n_OL = numel(velM);
%% CL, only create design matrix to use for prediction

velM_cl = preprocessSignal(dataCL.vr_data_resampled.velM',21);
velP_cl = preprocessSignal(dataCL.vr_data_resampled.velP',21);
pa = removeOutliers(dataCL.pupil_area,[1 99]);
pupil_area_cl = preprocessSignal(pa,11);

p_x = removeOutliers(dataCL.pupil_com(:,1),[1 99]);
pupil_x_cl = preprocessSignal(p_x,11);
p_y = removeOutliers(dataCL.pupil_com(:,2),[1 99]);
pupil_y_cl = preprocessSignal(p_y,11);

%%



% take position mod length of track (AFTER computing speed)


glmData=struct();
sample_idx = 1:downsample:numel(velM);
%% create A

var_name = {'velM','velP','pupilArea','pupilX','pupilY'};
A={};
s = 0.5; % spline parameter

%% velM and velP
spdVec = linspace(min([velM; velP;velM_cl;]),max([velM; velP; velM_cl]),15);
spdVec(1)=spdVec(1)-0.1;
spdVec(end)=spdVec(end)+0.1;
%spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[velMgrid,ctl_pts_velM] = spline_1d([velM;velM_cl],spdVec,s);
[velPgrid,ctl_pts_velP] = spline_1d([velP; velP_cl],spdVec,s);


A{1} = velMgrid(1:n_OL,:);
B{1} = velMgrid(n_OL+1:end,:);
all_control_points{1} = ctl_pts_velM;
A{2} = velPgrid(1:n_OL,:);
B{2} = velPgrid(n_OL+1:end,:);
all_control_points{2} = ctl_pts_velP;
%% pupil
areaVec = linspace(min([pupil_area;pupil_area_cl]),max([pupil_area;pupil_area_cl]),11);
areaVec(1)=[areaVec(1)-0.01];
areaVec(end)=areaVec(end)+0.01;
[areaGrid,ctl_pts_area] = spline_1d([pupil_area;pupil_area_cl],areaVec,s);

xVec = linspace(min([pupil_x;pupil_x_cl]),max([pupil_x;pupil_x_cl]),11);
xVec(1)=xVec(1)-0.01;
xVec(end)=xVec(end)+0.01;
[xGrid,ctl_pts_x] = spline_1d([pupil_x;pupil_x_cl],xVec,s);

yVec = linspace(min([pupil_y;pupil_y_cl]),max([pupil_y;pupil_y_cl]),11);
yVec(1)=yVec(1)-0.01;
yVec(end)=yVec(end)+0.01;
[yGrid,ctl_pts_y] = spline_1d([pupil_y;pupil_y_cl],yVec,s);

A{3}=areaGrid(1:n_OL,:);
B{3}=areaGrid(n_OL+1:end,:);
A{4}=xGrid(1:n_OL,:);
B{4}=xGrid(n_OL+1:end,:);
A{5}=yGrid(1:n_OL,:);
B{5}=yGrid(n_OL+1:end,:);
all_control_points{3} = ctl_pts_area;
all_control_points{4} = ctl_pts_x;
all_control_points{5} = ctl_pts_y;

%%
for iA=1:numel(A)
    A{iA}=A{iA}(sample_idx,:);
end
%% MM
%A{3}=bl_trial;

%%
for cellIDX=1:length(good_cells)
    spike_t = dataOL.sp.st(dataOL.sp.clu==good_cells(cellIDX));
    %[~,~,spike_idx] = histcounts(spike_t,post);
    
    [spiketrain,~,spikeIDX] = histcounts(spike_t,dataOL.post(sample_idx));
    spiketrain = spiketrain';
    numFolds = 10;
    T = numel(spiketrain);
    numPts = 10*round(1/opt.TimeBin); % 3 seconds. i've tried #'s from 1-10 seconds.. not sure what is best
    for iA=1:numel(A)
        A{iA}=A{iA}(1:T,:);
    end
    
    %%
    [train_ind,test_ind] = compute_test_train_ind(numFolds,numPts,T);
    %     trial = dataOL.trial(sample_idx);
    %     tmp_ind = mod(trial(1:T),numFolds);
    %     for k = 1:numFolds
    %         test_ind{k} = find(tmp_ind == k-1);
    %         train_ind{k} = setdiff(1:T,test_ind{k});
    %     end
    %%%%%%%% FORWARD SEARCH PROCEDURE %%%%%%%%%
    mean_fr = mean(spiketrain); %mean per bin
    yhat =mean_fr*ones(size(A{1},1),1);
        yhat_cl = mean_fr*ones(size(B{1},1),1);
    try
        fprintf('\t Fitting model  for cell %d \n', cellIDX);
        [allModelTestFits, allModelTrainFits, bestModels, bestModelFits, parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,train_ind,test_ind);
        
        vars = sort(bestModels);
        X=ones(length(spiketrain),1);
        BX = ones(length(velM_cl),1);
        for iV=vars
            X=[X A{iV}];
            BX = [BX B{iV}];
        end
        if final_pval<=0.05 %only use model if fit is sig better than average firing rate model
        yhat = exp(X*parameters{end}');
        yhat_cl = exp(BX*parameters{end}');
        end
        
    catch ME
        fprintf('Model fitting failed for %d \n',cellIDX)
        warning(ME.message)
        allModelTestFits = NaN;
        bestModels = NaN;
        bestModelFits = NaN;
        parameters = NaN;
        pvals = [];
        tuning_curves ={};
        final_pval = nan;
        
    end
    try
        tuning_curves = glm_tuning_curves_noPos(A,bestModels,parameters{end},all_control_points,s,median(diff(dataOL.post(sample_idx))));
    catch ME
        fprintf('Tuning curves failed \n')
    end
    
    
    
    
    
    glmData(cellIDX).allModelTestFits = allModelTestFits;
    glmData(cellIDX).bestModels = bestModels;
    glmData(cellIDX).bestModelFits = bestModelFits;
    glmData(cellIDX).parameters = parameters;
    glmData(cellIDX).pvals = pvals;
    glmData(cellIDX).yhat = yhat;
    glmData(cellIDX).yhat_cl = yhat_cl;
    glmData(cellIDX).var_name = var_name;
    glmData(cellIDX).final_pval = final_pval;
    glmData(cellIDX).tuning_curves  = tuning_curves;
    
end

function signal = preprocessSignal(signal,smoothWin)
signal = fillmissing(signal,'nearest');
fi = gausswin(smoothWin);
fi = fi/sum(fi);
signal = conv([repmat(signal(1),floor(smoothWin/2),1); signal; repmat(signal(end),floor(smoothWin/2),1)],fi,'valid');


function signal = removeOutliers(signal,quantiles)
quants=prctile(signal,quantiles);
mi = quants(1);
ma = quants(2);
signal(signal<mi)=mi;
signal(signal>ma)=ma;


