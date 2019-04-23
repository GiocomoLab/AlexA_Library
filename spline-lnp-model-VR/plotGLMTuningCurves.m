function [ fighandle ] = plotGLMTuningCurves(fighandle,var_name,posx,speed, max_speed,bl_posx,A,s,all_control_points,spiketrain,spike_idx,trial,params,bestModels,allModelTestFits,parameters)
    figure(fighandle)
    num_plot_columns=3;
    numVar=length(var_name);
    %spatial firing map
    subplot(4,3,[2 5 8]);
    plot(posx(spike_idx),trial(spike_idx),'k.');
    
    ylim([0 max(trial)+1]); xlim([params.TrackStart params.TrackEnd]);
    title(['best Models: ' num2str(bestModels)])
    
    xlabel('pos'); ylabel('trial');
    %posbins = linspace(params.TrackStart,params.TrackEnd,40);
    %[a,b,c]=histcounts(po
    [pos_tuning_curve,pos_occupancy,bins] = compute_1d_tuning_curve(bl_posx,spiketrain,30,params.TrackStart,params.TrackEnd);
    xbincent=0.5 * (bins(1:end-1) + bins(2:end));
    subplot(4,num_plot_columns,1);
    plot(xbincent,pos_occupancy.*params.TimeBin)
    title('position occupancy')
    axis tight
    
    subplot(4,num_plot_columns,1+num_plot_columns);
    plot(xbincent,pos_tuning_curve./params.TimeBin)
    title('position tuning curve')
    ylabel('spikes/s')
    axis tight
    
    [speed_tuning_curve,speed_occupancy,bins] = compute_1d_tuning_curve(speed,spiketrain,10,0,max_speed);
    
    subplot(4,num_plot_columns,3);
    plot(linspace(0,max_speed,10),speed_occupancy.*params.TimeBin,'k','linewidth',2)
    box off
    title('speed occupancy')
    axis tight
    ylabel('seconds')
    
    subplot(4,num_plot_columns,3+num_plot_columns);
    plot(.5*bins(1:end-1)+.5*bins(2:end),speed_tuning_curve./params.TimeBin,'k','linewidth',2)
    box off
    title('speed tuning curve')
    axis tight
    ylabel('spikes/s')
    
    
    plotfig = 1;
    final_param = parameters{end};
    [tuning_curves] = plot_all_tuning(A,bestModels,final_param,all_control_points,s,fighandle,params.TimeBin);
    
    
    firstModelFit = allModelTestFits{1};
    fig1 = subplot(4,num_plot_columns,num_plot_columns*3+1:num_plot_columns*3+num_plot_columns);
    errorbar(1:numVar,mean(firstModelFit),std(firstModelFit)/sqrt(10),'.k','linewidth',2)
    hold on
    plot([1 numVar],[0 0],'--b','linewidth',1);
    hold off
    box off
    set(gca,'xtick',1:length(var_name))
    set(gca,'xticklabel',var_name)
    ylabel('bits/spike')


end

