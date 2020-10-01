function [glmData] = fitGLM_mismatch(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here





   
    %bl_trial=data.trial(bl_idx);
    params.TimeBin = 0.02;
    speed=data.true_speed;
speed(speed<0) = 0;
max_speed=ceil(prctile(speed,99));
speed(speed > max_speed) = max_speed;
% take position mod length of track (AFTER computing speed)
posx = data.posx;
posx(posx<0)=0;
posx(posx>400)=400;

glmData=struct();
keep_data={};
%% create A

var_name = {'position','speed'};
A={};
all_control_points={};
s = 0.5; % spline parameter
%% position
numctrlpoints_pos=30;
x_vec = linspace(0,400,numctrlpoints_pos);
x_vec(1) = x_vec(1)-0.01;

[posgrid,ctl_pts_pos] = spline_1d(posx(2:end),x_vec,s);
A{1} = posgrid;
all_control_points{1} = ctl_pts_pos;

%% speed
spdVec = linspace(0,max_speed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[speedgrid,ctl_pts_speed] = spline_1d(speed(2:end),spdVec,s);
A{2} = speedgrid;
all_control_points{2} = ctl_pts_speed;
%% gain
%A{3}=bl_trial;

%%
good_cells = data.sp.cids(data.sp.cgs==2);
for cellIDX=17:length(good_cells)
    spike_t = data.sp.st(data.sp.clu==good_cells(cellIDX));
    %[~,~,spike_idx] = histcounts(spike_t,post);
    
    spiketrain = histcounts(spike_t,data.post)';
    numFolds = 10;
    T = numel(spiketrain);
    numPts = 3*round(1/params.TimeBin); % 3 seconds. i've tried #'s from 1-10 seconds.. not sure what is best
    
    
    %%
    [train_ind,test_ind] = compute_test_train_ind(numFolds,numPts,T);
    
    %%%%%%%% FORWARD SEARCH PROCEDURE %%%%%%%%%
    try 
        fprintf('\t Fitting model  for cell %d \n', cellIDX);
        [allModelTestFits, allModelTrainFits, bestModels, bestModelFits, parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,train_ind,test_ind);

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
        tuning_curves = glm_tuning_curves(A,bestModels,parameters{end},all_control_points,s,0.02);
    catch ME
        fprintf('Tuning curves failed \n')
    end

    glmData(cellIDX).allModelTestFits = allModelTestFits;
    glmData(cellIDX).bestModels = bestModels;
    glmData(cellIDX).bestModelFits = bestModelFits;
    glmData(cellIDX).parameters = parameters;
    glmData(cellIDX).pvals = pvals;
    glmData(cellIDX).all_control_points = all_control_points;
    glmData(cellIDX).var_name = var_name;
    glmData(cellIDX).final_pval = final_pval;
    glmData(cellIDX).tuning_curves  = tuning_curves;

end

