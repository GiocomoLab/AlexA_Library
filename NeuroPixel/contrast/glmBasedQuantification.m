%contrast = load('Z:\giocomo\attialex\NP_DATA\npJ4_0512_contrast_1.mat');
A={};
all_control_points={};
s = 0.5; % spline parameter

numctrlpoints_pos=30;
x_vec = linspace(0,400,numctrlpoints_pos);
x_vec(1) = x_vec(1)-0.01;
posx=contrast.posx;
posx(posx<0)=0;
posx(posx>400)=400;
[posgrid,ctl_pts_pos] = spline_1d(posx,x_vec,s);
all_control_points{1} = ctl_pts_pos;
good_cells = contrast.sp.cids(contrast.sp.cgs==2);

% evaluate model
error_ratio = nan(numel(good_cells),max(contrast.trial));
for cellID=1:length(good_cells)
    if  numel(contrast.glmData(cellID).bestModels)==1 && contrast.glmData(cellID).bestModels(1)==1 && mean(contrast.glmData(cellID).bestModelFits)>0
        parameters = contrast.glmData(cellID).parameters{end};
        b0 = parameters(1);
        param = parameters(2:end);
        
        param1 = param;
        scale(1) = mean(exp(posgrid*param1'));
        fr = exp(posgrid*param1')*exp(b0)/0.02;
        
        %double check
        [~,baseline_s]=get_spatial_map(contrast,find(contrast.trial_contrast==100));
        [~,all_s] = get_spatial_map(contrast);
        bins = [0:4:400];
        bins(1)=-0.01;
        posx=contrast.posx;
        posx(posx<0)=0;
        posx(posx>400)=400;
        discrete_pos = discretize(contrast.posx,bins);
        tuning_curve=zeros(size(bins));
        for ib=1:length(bins);tuning_curve(ib)=mean(fr(discrete_pos==ib));end
        tuning_curve=tuning_curve(1:end-1);
        figure
        plot(tuning_curve)
        tc=mean(baseline_s(cellID,:,:),3);
        hold on
        plot(tc*50)
        legend('LNP Prediction','Actual Data')
        
        % evaluate baseline model
        spike_idx = contrast.sp.clu==good_cells(cellID);
        [~,~,c]=histcounts(contrast.sp.st(spike_idx),contrast.post);
        %time in baseline
        time_in_baseline=nnz(ismember(contrast.trial,find(contrast.trial_contrast ==100)))/50;
        spikes_in_baseline = nnz(ismember(contrast.trial(c),find(contrast.trial_contrast==100)));
        baseline_fr = spikes_in_baseline/time_in_baseline;
        
        for trial=1:max(contrast.trial)
            trial_idx = contrast.trial==trial;
            [o,~,~]=histcounts(posx(trial_idx),bins);
            baseline_firing = o*0.02*baseline_fr;
            current_spatial_firing=squeeze(all_s(cellID,:,trial))*50;
            error_model = current_spatial_firing-tuning_curve;
            error_model = sum(sqrt(error_model.^2));
            error_baseline = current_spatial_firing-baseline_firing;
            error_baseline = sum(sqrt((error_baseline).^2));
            error_ratio(cellID,trial)=(error_baseline-error_model)/(error_baseline+error_model);
        end
    end
end
    % figure
    % for trial=1:50
    %     subplot(2,1,1)
    %     trial_idx = contrast.trial==trial;
    %     %plot model prediction
    %     plot(contrast.post(trial_idx),fr(trial_idx))
    %     %plot spikes
    %     [~,~,c]=histcounts(contrast.sp.st(spike_idx),contrast.post);
    %     hold on
    %     c_idx = contrast.trial(c)==trial;
    %     scatter(contrast.post(c(c_idx)),ones(1,nnz(c_idx)),2)
    %     %spikes=contrast.sp.st(contrast.sp.clu==good_cells(3))
    %
    %     [o,~,~]=histcounts(posx(trial_idx),bins);
    %     baseline_firing = o*0.02*baseline_fr;
    %
    %     current_spatial_firing=squeeze(baseline_s(cellID,:,trial))*50;
    %     error_model = current_spatial_firing-tuning_curve;
    %     error_model = sum(sqrt(error_model.^2));
    %     error_baseline = current_spatial_firing-baseline_firing;
    %     error_baseline = sum(sqrt((error_baseline).^2));
    %
    %     subplot(2,1,2)
    %     plot(current_spatial_firing)
    %     hold on
    %     plot(tuning_curve)
    %
    %     plot(baseline_firing)
    %
    %     if(error_model>error_baseline)
    %         title('baseline')
    %     else
    %         title('model')
    %     end
    %     pause
    %     clf
    % end
    %
    
    %%
    
    figure
    for trial=1:50
        subplot(2,1,1)
        trial_idx = contrast.trial==trial;
        %plot model prediction
        plot(contrast.post(trial_idx),fr(trial_idx))
        %plot spikes
        [~,~,c]=histcounts(contrast.sp.st(spike_idx),contrast.post);
        hold on
        c_idx = contrast.trial(c)==trial;
        scatter(contrast.post(c(c_idx)),ones(1,nnz(c_idx)),2)
        %spikes=contrast.sp.st(contrast.sp.clu==good_cells(3))
    
        [o,~,~]=histcounts(posx(trial_idx),bins);
        baseline_firing = o*0.02*baseline_fr;
    
        current_spatial_firing=squeeze(baseline_s(cellID,:,trial))*50;
        error_model = current_spatial_firing-tuning_curve;
        error_model = sum(sqrt(error_model.^2));
        error_baseline = current_spatial_firing-baseline_firing;
        error_baseline = sum(sqrt((error_baseline).^2));
    
        subplot(2,1,2)
        plot(current_spatial_firing)
        hold on
        plot(tuning_curve)
    
        plot(baseline_firing)
    
        if(error_model>error_baseline)
            title('baseline')
        else
            title('model')
        end
        pause
        clf
    end
    