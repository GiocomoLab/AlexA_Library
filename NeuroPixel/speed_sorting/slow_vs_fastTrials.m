function slow_vs_fastTrials(filepath)

[~,session_name,~]=fileparts(filepath);

image_save_dir = strcat('/oak/stanford/groups/giocomo/attialex/Images/','speed_sort4');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
run('/home/users/attialex/AlexA_Library/default_paths.m')



load(filepath);
addpath(genpath('/home/users/attialex/AlexA_Library'));
params = readtable('UniversalParams.xlsx');
if ~isvarname('anatomy')
    return
end
try

if isfield(anatomy,'parent_shifted')
    region = anatomy.parent_shifted;
else
    region = anatomy.cluster_parent;
end
catch
    disp('no anatomy')
    return
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
%spatialMap=spatialMap(:,1:end-1,:);
%dwell_time=dwell_time(:,1:end-1);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
%% calc speed differential
speed = calcSpeed(posx,params);
posWindow= [10 390];
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


bl_trials = find(trial_gain == 1 & trial_contrast == 100);


bl_sorted = sid;
bl_sorted(~ismember(bl_sorted,bl_trials))=[];
nT=numel(bl_sorted);
M=round(.1*nT);

slow_idx = bl_sorted(1:M);
fast_idx = bl_sorted(nT-M:nT);
%h=figure('visible','off');

gains_all = [0.8 0.7 0.6 0.5 0.2];
gains = sort(unique(trial_gain),'descend');
gains = gains(2:end);


[~,gain_plot_idx] = ismember(gains,gains_all);
plot_colors_gain = cool(numel(gains_all));
fig=figure('visible','on');
delay_per_cell = zeros(2,numel(good_cells));
trials_fast = zeros(numel(edges)-1,numel(good_cells));
trials_slow = trials_fast;
for cellIDX = 1:numel(good_cells)
    cluID = find(sp.cids==good_cells(cellIDX));
    clu_reg{cellIDX}=region{cluID};
    mS=spMap(:,:,cellIDX);
    mS(isnan(mS))=0;
    mS=conv2(mS,filt,'same');
    tmp_slow = nanmean(mS(:,slow_idx),2);
    tmp_fast = nanmean(mS(:,fast_idx),2);
    trials_slow(:,cellIDX)=tmp_slow;
    trials_fast(:,cellIDX)=tmp_fast;
    [rr,lags]=xcorr(tmp_slow(posBin(1):posBin(2)),tmp_fast(posBin(1):posBin(2)),5,'coeff');
    [~,tmp]=max(rr);
    delay_per_cell(1,cellIDX)=lags(tmp)*mean(diff(edges));
    delay_per_cell(2,cellIDX)=rr(tmp);
    %dd=finddelay(tmp_slow(posBin(1):posBin(2)),tmp_fast((posBin(1):posBin(2))));
    ax=subplot(1,1,1);
    p1=plot(edges(2:end),(tmp_slow));
    hold(ax,'on')
    p2=plot(edges(2:end),tmp_fast);
    title(round(delay_per_cell(1,cellIDX)))

%     for ig=1:numel(gains)
%         
%         idx = trial_gain == gains(ig);
%         plot(nanmean(mS(:,idx),2),'Color',plot_colors_gain(gain_plot_idx(ig),:));
%     end
    for ij=posWindow
        xline(ij,'r');
    end
    legend([p1 p2],{'slow','fast'})

%     ax2=subplot(2,1,2);
%     hold(ax2,'on')
%     % get spike times and index into post
%     spike_t = sp.st(sp.clu==good_cells(cellIDX));
%     [~,~,spike_idx] = histcounts(spike_t,post);
%     
%     
%     
%     % gain trials
%     plot(posx(spike_idx),trial_sorted(spike_idx),'.','Color',[0 0 0])
%     for j = 1:numel(gains)
%         keep = trial_gain(trial(spike_idx))==gains(j);
%         plot(posx(spike_idx(keep)),trial_sorted(spike_idx(keep)),'.','Color',plot_colors_gain(gain_plot_idx(j),:));
%     end
%     xlim([params.TrackStart params.TrackEnd]);
%     ylim([0 max(trial)+1]);
%     title(sprintf('c%d, %s',good_cells(cellIDX),region{cluID}));
%     xticks(''); yticks('');
%     for ij=[80 :80:320]
%         xline(ij);
%     end
%     for ij=posWindow
%         xline(ij,'r');
%     end
    
    
    %saveas(fig,fullfile(image_save_dir,sprintf('%s_%s_%d.png',clu_reg{cellIDX},session_name,cellIDX)),'png');
    clf
    %%
end

data.delay = delay_per_cell;

data.region = clu_reg;
data.session = session_name;
data.slow = mean(trial_speed(slow_idx));

data.fast = mean(trial_speed(fast_idx));
data.slow_trials = trials_slow;
data.fast_trials = trials_fast;

save(fullfile(OAK,'attialex','speed_sort4',session_name),'data')