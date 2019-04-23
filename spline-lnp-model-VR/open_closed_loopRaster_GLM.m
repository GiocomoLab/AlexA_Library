load('Z:\giocomo\attialex\NP_DATA\mismatch\npH3_0404_playback_1.mat')
open_loop=load('Z:\giocomo\attialex\NP_DATA\mismatch\npH3_0404_playback_1.mat')
closed_loop=load('Z:\giocomo\attialex\NP_DATA\mismatch\npH3_0404_mismatch_1.mat')
%%
session_name='npH3_0404_playback_1';
image_save_dir = strcat('F:\Data\npAna\images\',session_name,'\pretty_rasters\');
if ~isdir(image_save_dir)
    mkdir(image_save_dir);
end

%%
run_speed = gauss_smoothing(true_speed/params.TimeBin,10);
vis_speed = calcSpeed(posx,params);
run_speed(run_speed<0) = 0;
vis_speed(vis_speed<0) = 0;
max_runspeed=ceil(prctile(run_speed,99));
max_visspeed=ceil(prctile(vis_speed,99));
vis_speed(vis_speed > max_visspeed) = max_visspeed;
run_speed(run_speed > max_runspeed) = max_runspeed;
% take position mod length of track (AFTER computing speed)
posx(posx<0)=0;
posx(posx>params.TrackEnd)=params.TrackEnd;

%% create A

var_name = {'position','run_speed','vis_speed'};
A={};
all_control_points={};
s = 0.5; % spline parameter

numctrlpoints_pos=30;
x_vec = linspace(params.TrackStart,params.TrackEnd,numctrlpoints_pos);
x_vec(1) = x_vec(1)-0.01;
[posgrid,ctl_pts_pos] = spline_1d(posx,x_vec,s);
A{1} = posgrid;
all_control_points{1} = ctl_pts_pos;


spdVec = linspace(0,max_runspeed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[speedgrid,ctl_pts_speed] = spline_1d(run_speed,spdVec,s);
A{2} = speedgrid;
all_control_points{2} = ctl_pts_speed;


spdVec = linspace(0,max_visspeed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[vis_speedgrid,ctl_pts_visspeed] = spline_1d(vis_speed,spdVec,s);
A{3} = vis_speedgrid;
all_control_points{3} = ctl_pts_visspeed;
%%
good_cells = sp.cids(sp.cgs==2);
        fighandle = figure('Position',[680   189   728   789]);

for cellIDX=1:length(good_cells)

    spike_t = sp.st(sp.clu==good_cells(cellIDX));
    [~,~,spike_idx] = histcounts(spike_t,post);
    ax=subplot(4,2,[1 3 5 7]);
    hold(ax);
    spatial_raster(closed_loop.sp,good_cells(cellIDX),closed_loop.post,closed_loop.posx,closed_loop.trial,ax,'b');
    offset=max(closed_loop.trial);
    spatial_raster(open_loop.sp,good_cells(cellIDX),open_loop.post,open_loop.posx,open_loop.trial+offset,ax,'r');
    spiketrain = histc(spike_t,post);
    
    xlim([params.TrackStart params.TrackEnd])
    
    plotTuningCurves(subplot(4,2,2),subplot(4,2,2),spiketrain,'position',posx,30,0,400,params.TimeBin)
    plotTuningCurves(subplot(4,2,4),subplot(4,2,4),spiketrain,'run speed',run_speed,10,0,max_runspeed,params.TimeBin)
    plotTuningCurves(subplot(4,2,6),subplot(4,2,6),spiketrain,'vis speed',vis_speed,10,0,max_visspeed,params.TimeBin)

    for iV=closed_loop.glmData(cellIDX).bestModels
        subplot(4,2,iV*2)
        hold on
        plot(closed_loop.glmData(cellIDX).glm_xVals{iV},closed_loop.glmData(cellIDX).glmTuningCure{iV},'b--')
    end
    
    for iV=open_loop.glmData(cellIDX).bestModels
        subplot(4,2,iV*2)
        hold on
        plot(open_loop.glmData(cellIDX).glm_xVals{iV},open_loop.glmData(cellIDX).glmTuningCure{iV},'r--')
    end
    
    
    firstModelFit = closed_loop.glmData(cellIDX).allModelTestFits{1};
    subplot(4,2,8)
    errorbar(1:numel(var_name),mean(firstModelFit),std(firstModelFit)/sqrt(10),'.b','linewidth',2)
    hold on
    plot([1 numel(var_name)],[0 0],'--b','linewidth',1);
    box off
    set(gca,'xtick',1:length(var_name))
    set(gca,'xticklabel',var_name,'TickLabelInterpreter','None')
    ylabel('bits/spike')
    
    firstModelFit = open_loop.glmData(cellIDX).allModelTestFits{1};
    errorbar(1:numel(var_name),mean(firstModelFit),std(firstModelFit)/sqrt(10),'.r','linewidth',2)
    
    saveas(fighandle,fullfile(image_save_dir,sprintf('%d.png',good_cells(cellIDX))),'png');

    clf
end