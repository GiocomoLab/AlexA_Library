% plots rasters for all cells in all sessions
% colors spikes by gain and contrast
% MGC 2/28/2019

%% params
% make sure paths are correct
addpath(genpath('/home/users/attialex/AlexA_Library'))
% some params
params = readtable('UniversalParams.xlsx');
save_figs = true;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
%data_dir = 'Z:\giocomo\attialex\NP_DATA\';
data_dir = '/oak/stanford/groups/giocomo/attialex/NP_DATA';

%session_name = {'AA5_190809_gain_1'};
session_name = {};
sn = dir(fullfile(data_dir,'AA*.mat'));
for iS = 1:numel(sn)
    if contains(sn(iS).name,'gain')  %~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback'))
        session_name{end+1}=sn(iS).name(1:end-4);
    end
end
% all the values we ever use
gains_all = [0.8 0.7 0.6 0.5 0.2];
contrasts_all = [100 50 20 10 5 2 0];

%% iterate over sessions
for session_num = 1:numel(session_name)    
    % load data
    fprintf('session %d/%d: %s\n',session_num,numel(session_name),session_name{session_num});
    clear anatomy
    load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));
    
    cells_to_plot = sp.cids(sp.cgs==2); % all good cells
    
    % make image dir if it doesn't exist
    %image_save_dir = fullfile('F:','images','rastersRegion');
    image_save_dir = fullfile('/oak/stanford/groups/giocomo/attialex/images/rastersRegionParent_subset20');
    if ~isfolder(image_save_dir)
        mkdir(image_save_dir)
    end
    
    if ~exist('anatomy')
        continue
    end
    
    if isfield(anatomy,'region_shifted')
        %regionID=anatomy.region_shifted;
        regionID = anatomy.parent_shifted;
    else
        %regionID = anatomy.cluster_region;
        regionID = anatomy.cluster_parent;
    end
    regionID = regionID(sp.cgs==2);

    % get plot colors for gains
    gains = sort(unique(trial_gain),'descend');
    contrasts = sort(unique(trial_contrast),'descend');
    gains = gains(2:end);
    if contains(session_name{session_num},'dark')
       gains = [];
       contrasts = 100;
       trial_gain = 1*ones(1,max(trial));
       trial_contrast = 100*ones(1,max(trial));
    end
    
    [~,gain_plot_idx] = ismember(gains,gains_all);
    plot_colors_gain = cool(numel(gains_all));

    % get plot colors for contrasts
    
    [~,contrast_plot_idx] = ismember(contrasts,contrasts_all);
    plot_colors_contrast = gray(numel(contrasts_all)+1);
    plot_colors_contrast = plot_colors_contrast(1:end-1,:);


    spike_depth = anatomy.tip_distance;
    
    %% make raster plots for all cells
    set(0,'DefaultFigureRenderer','painters')
    h = figure('Position',[100   467   224   133]); hold on;
    
    for k = 1:numel(cells_to_plot)
        
        cla;
        fprintf('\tcell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

        % get spike times and index into post
        spike_t = sp.st(sp.clu==cells_to_plot(k));
        [~,~,spike_idx] = histcounts(spike_t,post);

        % make figure


        % plot spike raster
        % baseline trials for all contrasts
        for j = 1:numel(contrasts)
            keep = trial_contrast(trial(spike_idx))==contrasts(j) & ...
                trial_gain(trial(spike_idx))==1;
            plot(posx(spike_idx(keep)),trial(spike_idx(keep)),'.','Color',[0 0 0]);
        end
        % gain trials
        for j = 1:numel(gains)
            keep = trial_gain(trial(spike_idx))==gains(j);
            plot(posx(spike_idx(keep)),trial(spike_idx(keep)),'.','Color',plot_colors_gain(gain_plot_idx(j),:));
        end
        xlim([params.TrackStart params.TrackEnd]);
        ylim([0 max(trial)+1]);
        %title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
        xticks(''); yticks('');
        %
        rID = regionID{k};
        rID = strrep(rID,'/','_');
        rID = strrep(rID,'\','_');
        if ~isfolder(fullfile(image_save_dir,rID))
            mkdir(fullfile(image_save_dir,rID));
        end
        % save fig
        set(gca,'YLim',[0 34])
        %axis ij
        if save_figs
            saveas(h,fullfile(image_save_dir,rID,sprintf('%s_%d.pdf',session_name{session_num},k)),'pdf');
            %print(fullfile(image_save_dir,rID,sprintf('%s_%d.png',session_name{session_num},k)),'-dpng','-r600')
            %saveas(h,fullfile(image_save_dir,rID,sprintf('%d.pdf',k)),'pdf');
        end
        cla;
%         catch ME
%             warning(sprintf('did not work for cell %d',k))
    
    %end
    
        
    end
end