%% preprocess variables
speed = calcSpeed(posx,params);
speed(speed<0) = 0;
max_speed=ceil(prctile(speed,99));
speed(speed > max_speed) = max_speed;
% take position mod length of track (AFTER computing speed)
posx(posx<0)=0;
posx(posx>params.TrackEnd)=params.TrackEnd;
%% create A

var_name = {'position','speed'};
A={};
all_control_points={};
s = 0.5; % spline parameter
%% position
numctrlpoints_pos=30;
x_vec = linspace(params.TrackStart,params.TrackEnd,numctrlpoints_pos);
x_vec(1) = x_vec(1)-0.01;
[posgrid,ctl_pts_pos] = spline_1d(posx,x_vec,s);
A{1} = posgrid;
all_control_points{1} = ctl_pts_pos;

%% speed
spdVec = linspace(0,max_speed,15);
spdVec(1)=-0.01; %for boundary condition
%maybe use linspace?
[speedgrid,ctl_pts_speed] = spline_1d(speed,spdVec,s);
A{2} = speedgrid;
all_control_points{2} = ctl_pts_speed;

%%
good_cells = sp.cids(sp.cgs==2);
spike_t = sp.st(sp.clu==good_cells(2));
[~,~,spike_idx] = histcounts(spike_t,post);

spiketrain = histc(spike_t,post);
numFolds = 10;
T = numel(spiketrain);
numPts = 3*round(1/params.TimeBin); % 3 seconds. i've tried #'s from 1-10 seconds.. not sure what is best


%%
[train_ind,test_ind] = compute_test_train_ind(numFolds,numPts,T);

%%%%%%%% FORWARD SEARCH PROCEDURE %%%%%%%%%
[allModelTestFits, allModelTrainFits, bestModels, bestModelFits, parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,train_ind,test_ind);

%% some plotting
num_plot_columns=3;
numVar=length(var_name);
figure
%spatial firing map
subplot(4,3,[2 5 8]);
plot(posx(spike_idx),trial(spike_idx),'k.');
ylim([0 max(trial)+1]); xlim([params.TrackStart params.TrackEnd]);
xlabel('pos'); ylabel('trial');
%posbins = linspace(params.TrackStart,params.TrackEnd,40);
%[a,b,c]=histcounts(po
[pos_tuning_curve,pos_occupancy,bins] = compute_1d_tuning_curve(posx,spiketrain,30,params.TrackStart,params.TrackEnd);
xbincent=0.5 * (bins(1:end-1) + bins(2:end));
fig1 = subplot(4,num_plot_columns,1);
plot(xbincent,pos_occupancy.*params.TimeBin)
title('position occupancy')
axis tight

fig1 = subplot(4,num_plot_columns,1+num_plot_columns);
plot(xbincent,pos_tuning_curve./params.TimeBin)
title('position tuning curve')
ylabel('spikes/s')
axis tight

[speed_tuning_curve,speed_occupancy,bins] = compute_1d_tuning_curve(speed,spiketrain,10,0,max_speed);

fig1 = subplot(4,num_plot_columns,3);
plot(linspace(0,max_speed,10),speed_occupancy.*params.TimeBin,'k','linewidth',2)
box off
title('speed occupancy')
axis tight
ylabel('seconds')

fig1 = subplot(4,num_plot_columns,3+num_plot_columns);
plot(.5*bins(1:end-1)+.5*bins(2:end),speed_tuning_curve./params.TimeBin,'k','linewidth',2)
box off
title('speed tuning curve')
axis tight
ylabel('spikes/s')


plotfig = 1;
final_param = parameters{end};
[tuning_curves,fig1] = plot_all_tuning(A,bestModels,final_param,all_control_points,s,plotfig,params.TimeBin);


firstModelFit = allModelTestFits{1};
fig1 = subplot(4,num_plot_columns,num_plot_columns*3+1:num_plot_columns*3+num_plot_columns);
errorbar(1:numVar,mean(firstModelFit),std(firstModelFit)/sqrt(10),'.k','linewidth',2)
hold on
plot([1 numVar],[0 0],'--b','linewidth',1);
hold off
box off
set(gca,'xtick',[1 2])
set(gca,'xticklabel',var_name)
ylabel('bits/spike')
axis([0.5  2.5 -inf inf])

% save fig and close
    %saveas(fig1, sprintf('images/glm/%s/%d.png',session_name,good_cells(i)), 'png'); 
    %close(1);
