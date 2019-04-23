
    bl_posx=posx;
    bl_post=post;
    bl_spike_idx = true(size(sp.st));
    bl_trial = trial;
%TODO: adjust sp.st, sp.clu


%% preprocess variables
run_speed = gauss_smoothing(true_speed/params.TimeBin,10);
vis_speed = calcSpeed(posx,params);
run_speed(run_speed<0) = 0;
vis_speed(vis_speed<0) = 0;
max_runspeed=ceil(prctile(run_speed,99));
max_visspeed=ceil(prctile(vis_speed,99));
vis_speed(vis_speed > max_visspeed) = max_visspeed;
run_speed(run_speed > max_runspeed) = max_runspeed;
% take position mod length of track (AFTER computing speed)
bl_posx(bl_posx<0)=0;
bl_posx(bl_posx>params.TrackEnd)=params.TrackEnd;

glmData=struct();
keep_data={};
%% create A

var_name = {'position','run_speed','vis_speed'};
A={};
all_control_points={};
s = 0.5; % spline parameter
%% position
numctrlpoints_pos=30;
x_vec = linspace(params.TrackStart,params.TrackEnd,numctrlpoints_pos);
x_vec(1) = x_vec(1)-0.01;
[posgrid,ctl_pts_pos] = spline_1d(bl_posx,x_vec,s);
A{1} = posgrid;
all_control_points{1} = ctl_pts_pos;

%% run speed
spdVec = linspace(0,max_runspeed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[speedgrid,ctl_pts_speed] = spline_1d(run_speed,spdVec,s);
A{2} = speedgrid;
all_control_points{2} = ctl_pts_speed;
%%
%% run speed
spdVec = linspace(0,max_visspeed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[vis_speedgrid,ctl_pts_visspeed] = spline_1d(vis_speed,spdVec,s);
A{3} = vis_speedgrid;
all_control_points{3} = ctl_pts_visspeed;
%% gain
%A{3}=bl_trial;

%% 4
good_cells = sp.cids(sp.cgs==2);
for cellIDX=1:length(good_cells)
    spike_t = sp.st(sp.clu==good_cells(cellIDX) & bl_spike_idx);
    [~,~,spike_idx] = histcounts(spike_t,post);
    
    spiketrain = histc(spike_t,bl_post);
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
    end
    glmData(cellIDX).allModelTestFits = allModelTestFits;
    glmData(cellIDX).bestModels = bestModels;
    glmData(cellIDX).bestModelFits = bestModelFits;
    glmData(cellIDX).parameters = parameters;
    glmData(cellIDX).pvals = pvals;
    glmData(cellIDX).all_control_points = all_control_points;
    glmData(cellIDX).var_name = var_name;
    if numel(bestModels)>0
    [glmTuningCurve,glm_xVals]=getGLMTuningCurve(A,bestModels,parameters{end},all_control_points,s,params.TimeBin);
    glmData(cellIDX).glmTuningCure = glmTuningCurve;
    glmData(cellIDX).glm_xVals = glm_xVals;
    end
    
    %% some plotting
    if plotfig
    fighandle = figure('Position',[680   441   911   537]);
    plotTuningCurves(subplot(3,3,1),subplot(3,3,4),spiketrain,'position',posx,30,0,400,params.TimeBin)
    plotTuningCurves(subplot(3,3,2),subplot(3,3,5),spiketrain,'run speed',run_speed,10,0,max_runspeed,params.TimeBin)
    plotTuningCurves(subplot(3,3,3),subplot(3,3,6),spiketrain,'vis speed',vis_speed,10,0,max_visspeed,params.TimeBin)
    
    for iV=bestModels
        subplot(3,3,3+iV)
        hold on
        plot(glm_xVals{iV},glmTuingCurve_c{iV},'r--')
    end
    
    firstModelFit = allModelTestFits{1};
    subplot(3,3,7)
    errorbar(1:numel(var_name),mean(firstModelFit),std(firstModelFit)/sqrt(10),'.k','linewidth',2)
    hold on
    plot([1 numel(var_name)],[0 0],'--b','linewidth',1);
    hold off
    box off
    set(gca,'xtick',1:length(var_name))
    set(gca,'xticklabel',var_name,'TickLabelInterpreter','None')
    ylabel('bits/spike')
    
    saveas(fighandle, fullfile(plot_path,sprintf('glm_closedloop_%d.png',good_cells(cellIDX))));
    close(fighandle);
    end
end