%dark = load('Z:\giocomo\attialex\NP_DATA\npI1_0417_dark_1.mat')
%%
session_name='npH3_0401_dark_3';
image_save_dir = strcat('Z:\giocomo\attialex\images\',session_name,'\glm_dark\');
if ~isdir(image_save_dir)
    mkdir(image_save_dir);
end

%%
dp2=diff(dark.posx);
    teleport = find(dp2<-1);
    dp2(teleport)=0.5*dp2(teleport-1)+0.5*dp2(teleport+1);
    total_distance = [0; cumsum(dp2)];
    
    if min(diff(total_distance))<-1
        %seems some teleports happen across two frames, go again
        dp2=diff(dark,posx);
        teleport = find(dp2<-1);
        double_teleports = teleport((diff(teleport)==1));
        teleport(find(diff(teleport)==1)+1)=[];
        teleport(ismember(teleport,double_teleports))=[];
        dp2(teleport)=0.5*dp2(teleport-1)+0.5*dp2(teleport+1);
        for dp=1:length(double_teleports)
            tmp_i=double_teleports(dp);
            dp2([tmp_i, tmp_i+1])= 0.5*dp2(tmp_i-1)+0.5*dp2(tmp_i+2);
        end
        total_distance = [0; cumsum(dp2)];
        %single_teleports = teleport(diff(teleport)>1);
    end
    figure;plot(total_distance)
%%
params = readtable('UniversalParams.xlsx');
run_speed = gauss_smoothing([diff(total_distance)]/params.TimeBin,10);

run_speed(run_speed<0) = 0;
max_runspeed=ceil(prctile(run_speed,99));
run_speed(run_speed > max_runspeed) = max_runspeed;

run_bi = run_speed>1;
transitions=find(diff([0;run_bi;0]));
run_speed=[run_speed(1); run_speed];
run_times=transitions(2:2:end)-transitions(1:2:end);
%run_times(run_times<min_run)=[];
onsets=transitions(1:2:end);
offsets = transitions(2:2:end);

sit_time_transitions = find( diff([0;run_speed<1;0]));
sit_times = sit_time_transitions(2:2:end)-sit_time_transitions(1:2:end);
sit_offsets = sit_time_transitions(2:2:end);
[~,sit_sort]=sort(sit_times);

[a,b_sort]=sort(run_times);

run_times_sorted = run_times(b_sort);
onsets_sorted = onsets(b_sort);
offsets_sorted = offsets(b_sort);

run_numbers = zeros(size(total_distance));
time_since_onset = run_numbers;
distance_since_onset = run_numbers;

for iR=1:numel(run_times)
    current_onset = onsets_sorted(iR);
    current_length = run_times_sorted(iR);
    

        current_idx = current_onset:current_onset + current_length;
        run_numbers(current_idx) = iR;
        time_since_onset(current_idx) = dark.post(current_idx)-dark.post(current_idx(1));
        distance_since_onset(current_idx)=total_distance(current_idx)-total_distance(current_idx(1));

end
%%




%% create A

var_name = {'position','run_speed','distance_since_onset'};
A={};
all_control_points={};
s = 0.5; % spline parameter

numctrlpoints_pos=30;
x_vec = linspace(params.TrackStart,params.TrackEnd,numctrlpoints_pos);
x_vec(1) = x_vec(1)-0.01;
posx(posx>400)=400;
[posgrid,ctl_pts_pos] = spline_1d(posx,x_vec,s);
A{1} = posgrid;
all_control_points{1} = ctl_pts_pos;


spdVec = linspace(0,max_runspeed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[speedgrid,ctl_pts_speed] = spline_1d(run_speed,spdVec,s);
A{2} = speedgrid;
all_control_points{2} = ctl_pts_speed;

% create a for distance and time since onset
numctrlpoints_distance=30;
distance_x_vec = linspace(0,max(distance_since_onset),numctrlpoints_distance);
distance_since_onset(distance_since_onset<0)=0;
distance_x_vec(1)=-.01;

[distancegrid,ctl_pts_distance] = spline_1d(distance_since_onset,distance_x_vec,s);
A{3} = distancegrid;
all_control_points{3} = ctl_pts_distance;

%%
good_cells = sp.cids(sp.cgs==2);
        fighandle = figure('Position',[680   189   728   789]);
glmData = struct();
for cellIDX=44:length(good_cells)

    spike_t = sp.st(sp.clu==good_cells(cellIDX));
    [~,~,spike_idx] = histcounts(spike_t,post);
    ax=subplot(4,2,[1 3 5 7]);
    hold(ax);
    spatial_raster(dark.sp,good_cells(cellIDX),dark.post,dark.posx,dark.trial,ax,'b');
    
    spiketrain = histc(spike_t,post);
    numFolds = 10;
    T = numel(spiketrain);
    numPts = 10*round(1/params.TimeBin); % 3 seconds. i've tried #'s from 1-10 seconds.. not sure what is best
    
    
    
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
    if numel(bestModels)>0 && ~any(isnan(bestModels))
    [glmTuningCurve,glm_xVals]=getGLMTuningCurve(A,bestModels,parameters{end},all_control_points,s,params.TimeBin);
    glmData(cellIDX).glmTuningCure = glmTuningCurve;
    glmData(cellIDX).glm_xVals = glm_xVals;
    
    
    plotTuningCurves(subplot(4,2,2),subplot(4,2,2),spiketrain,'position',posx,30,0,400,params.TimeBin)
    plotTuningCurves(subplot(4,2,4),subplot(4,2,4),spiketrain,'run speed',run_speed,10,0,max_runspeed,params.TimeBin)
    plotTuningCurves(subplot(4,2,6),subplot(4,2,6),spiketrain,'distance',distance_since_onset,30,0,max(distance_since_onset),params.TimeBin)

    for iV=glmData(cellIDX).bestModels
        subplot(4,2,iV*2)
        hold on
        plot(glmData(cellIDX).glm_xVals{iV},glmData(cellIDX).glmTuningCure{iV},'b--')
    end
    
   
    
    
    firstModelFit = glmData(cellIDX).allModelTestFits{1};
    subplot(4,2,8)
    errorbar(1:numel(var_name),mean(firstModelFit),std(firstModelFit)/sqrt(10),'.b','linewidth',2)
    hold on
    plot([1 numel(var_name)],[0 0],'--b','linewidth',1);
    box off
    set(gca,'xtick',1:length(var_name))
    set(gca,'xticklabel',var_name,'TickLabelInterpreter','None')
    ylabel('bits/spike')
    
    firstModelFit = glmData(cellIDX).allModelTestFits{1};
    errorbar(1:numel(var_name),mean(firstModelFit),std(firstModelFit)/sqrt(10),'.r','linewidth',2)
    
    saveas(fighandle,fullfile(image_save_dir,sprintf('%d.png',good_cells(cellIDX))),'png');

    clf
    end
end