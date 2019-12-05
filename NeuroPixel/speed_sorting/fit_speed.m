function fit_speed(filepath)

[~,session_name,~]=fileparts(filepath);

plot_data=false;

run('/home/users/attialex/AlexA_Library/default_paths.m')

load(filepath);
addpath(genpath('/home/users/attialex/AlexA_Library'));
params = readtable('UniversalParams.xlsx');

if isfield(anatomy,'parent_shifted')
    region = anatomy.parent_shifted;
else
    region = anatomy.cluster_parent;
end

trials=[1:max(trial)];
%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
edges=[0:2:400];
edges(1)=-.01;
posx(posx<0)=0;
posx(posx>400)=400;
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.clu)+1);
    spatialMap=cat(3,spatialMap,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
good_cells = sp.cids(sp.cgs==2);

spatialMap=spatialMap(good_cells+1,:,:);
spatialMap=spatialMap(:,1:end-1,:);
dwell_time=dwell_time(:,1:end-1);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
%% calc speed differential
speed = calcSpeed(posx,params);
posWindow= [80 160];
%posWindow= [120 200];
posBin = [find(edges == posWindow(1)), find(edges == posWindow(2))];
trial_speed = zeros(1,max(trial));
for ii = 1:max(trial)
    idx = posx>posWindow(1) & posx<posWindow(2) & trial == ii;
    tmp = mean(speed(idx));
    
    trial_speed(ii)=tmp;
end
[~,sid]=sort(trial_speed);

trial_sorted = trial;
trial_rank = 1:max(trial);
trial_rank(sid)=trial_rank;
for ii=1:max(trial)
    idx = trial==ii;
    trial_sorted(idx)=trial_rank(ii);
end

spMap=shiftdim(spatialMap,1);



filt = gausswin(5);
filt = filt/sum(filt);

clu_reg = cell(numel(good_cells));
c_coeff = zeros(numel(good_cells),1);
p_vals = c_coeff;
models = cell(numel(good_cells));
if plot_data
    image_save_dir = fullfile('/oak/stanford/groups/giocomo/attialex/images/',...
        session_name,'speed_sort');
    if ~isfolder(image_save_dir)
        mkdir(image_save_dir)
    end
    
    
    fig = figure('Position',[109   487   742   259]); hold on;
    gains_all = [0.8 0.7 0.6 0.5 0.2];
    contrasts_all = [100 50 20 10 5 2 0];
    gains = sort(unique(trial_gain),'descend');
    contrasts = sort(unique(trial_contrast),'descend');
    gains = gains(2:end);
    
    
    [~,gain_plot_idx] = ismember(gains,gains_all);
    plot_colors_gain = cool(numel(gains_all));
    
    % get plot colors for contrasts
    
    [~,contrast_plot_idx] = ismember(contrasts,contrasts_all);
    plot_colors_contrast = gray(numel(contrasts_all)+1);
    plot_colors_contrast = plot_colors_contrast(1:end-1,:);
end

for cellIDX = 1:numel(good_cells)
    fprintf('now working on %d of %d \n',cellIDX,numel(good_cells))
    cluID = find(sp.cids==good_cells(cellIDX));
    clu_reg{cellIDX}=region{cluID};
    mS=spMap(:,:,cellIDX);
    mS(isnan(mS))=0;
    mS=conv2(mS,filt);
    
    [rr,lags]=xcorr(mS(posBin(1):posBin(2),:));
    nR=size(mS,2);
    delay_list=zeros(nR*(nR-1)/2,6); %delay, speed_diff, %i gain, i contrast, jgain j contrast
    idx = 0;
    nT = max(trial);
    for iT = 1:max(trial)
        for jT=(iT+1):max(trial)
            idx = idx+1;
            c_idx = (iT-1)*nT+jT;
            [mm,midx]=max(rr(:,c_idx));
            if ismember(midx,[1 size(rr,1)])
                midx = nan;
            end
            
            
            tmp_d = min(lags)-1;
            delay_list(idx,:)=[midx+tmp_d,trial_speed(iT)-trial_speed(jT), trial_gain(iT),trial_contrast(iT),trial_gain(jT),trial_contrast(jT) ];
        end
    end
    
    
    bl_trials = delay_list(:,3)==1 & delay_list(:,5)==1 & delay_list(:,4)==100 & delay_list(:,6)==100;
    mod = fitlm(delay_list(bl_trials,2),delay_list(bl_trials,1));
    [c,p]=corrcoef(delay_list(bl_trials,[2 1]),'Rows','complete');
    
    models{cellIDX}=mod.Coefficients;
    c_coeff(cellIDX)=c(2,1);
    p_vals(cellIDX)=p(2,1);
    
    if plot_data
        subplot(1,4,[1 2 3])
        hold on
        % get spike times and index into post
        spike_t = sp.st(sp.clu==good_cells(cellIDX));
        [~,~,spike_idx] = histcounts(spike_t,post);
        
        
        for j = 1:numel(contrasts)
            keep = trial_contrast(trial(spike_idx))==contrasts(j) & ...
                trial_gain(trial(spike_idx))==1;
            plot(posx(spike_idx(keep)),trial_sorted(spike_idx(keep)),'.','Color',plot_colors_contrast(contrast_plot_idx(j),:));
        end
        % gain trials
        for j = 1:numel(gains)
            keep = trial_gain(trial(spike_idx))==gains(j);
            plot(posx(spike_idx(keep)),trial_sorted(spike_idx(keep)),'.','Color',plot_colors_gain(gain_plot_idx(j),:));
        end
        xlim([params.TrackStart params.TrackEnd]);
        ylim([0 max(trial)+1]);
        title(sprintf('c%d, %s',good_cells(cellIDX),region{cluID}));
        xticks(''); yticks('');
        for ij=[80 :80:320]
            xline(ij);
        end
        for ij=posWindow
            xline(ij,'r');
        end
        
        
        
        subplot(1,4,4)
        plot(delay_list(bl_trials,2),delay_list(bl_trials,1),'.')
        axis square
        xlabel('speed diff')
        ylabel('delay')
        x1=min(delay_list(:,2));
        x2=max(delay_list(:,2));
        y1=[1 x1]*mod.Coefficients.Estimate;
        y2=[1 x2]*mod.Coefficients.Estimate;
        hold on
        plot([x1 x2],[y1 y2])
        
        saveas(fig,fullfile(image_save_dir,sprintf('%d.png',cellIDX)),'png');
        
        clf
    end
    
end
data.models = models;
data.c_coeff=c_coeff;
data.p_vals = p_vals;
data.region = clu_reg;
data.session = session_name;

save(fullfile(OAK,'attialex','speed_sort',session_name),'data')