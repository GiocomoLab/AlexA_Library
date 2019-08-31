%%
addpath(genpath('/home/users/attialex/AlexA_Library'));
addpath(genpath('/home/users/attialex/spikes/'));
%% find sessions based on table

session_info=readtable('/oak/stanford/groups/giocomo/attialex/NP_DATA/session_summary.xlsx');
nums = session_info(:,4);
nums=nums.Variables;
first_session = strcmp(nums,'1');
session_names = session_info(first_session,5);
session_names = session_names.Variables;
Files = struct();
for iF=1:numel(session_names)
    Files(iF).name=strcat(session_names{iF},'.mat');
end

%%
root='/oak/stanford/groups/giocomo/attialex/NP_DATA';
%root = 'Y:\giocomo\attialex\NP_DATA';
%root = fullfile('/oak/stanford/groups/giocomo','export','data','Projects','JohnKei_NPH3','E1','E1_190614_johncontrasttrack_train1_g0');
%Files = dir(fullfile(root,'*_contrast_*.mat'));
%Files = dir(fullfile(root,'*train*.mat'));

MERGED=struct;
cntr=1;
for iF=1:numel(Files)
    clearvars -except root Files iF MERGED cntr
    fn = fullfile(root,Files(iF).name);
    dataset = load(fn);
    try
    quantification_session
    MERGED(cntr).average_triggered = aa_spikes;
    MERGED(cntr).max_depth = maxChan_spikes;
    MERGED(cntr).name = Files(iF).name;
    MERGED(cntr).sp = dataset.sp; 
    cntr=cntr+1;
%     if mod(iF,7) ~=0
%         close(spikefig)
%     end
    catch ME
        ME.message
        warning(strcat('something wrong with ',Files(iF).name))
    end

end
%save(fullfile(root,'MERGED_FIRSTSESSIONS_DandA'),'MERGED');
%% scatter of max locations
bins = 0:40:3840;
time_bins = 0.002; 

tvec_spikes = [-100:100]*time_bins;

figure
hold on
col = summer(numel(MERGED));
for iF=1:numel(MERGED)
    [~,max_loc]=max(MERGED(iF).average_triggered,[],2);
    tmp_d=bins;
    tmp_d=tmp_d-bins(MERGED(iF).max_depth);
    plot(tvec_spikes(max_loc),tmp_d,'.','Color',col(iF,:))
        p=polyfit(tmp_d,tvec_spikes(max_loc),1);
    delay=polyval(p,tmp_d);
    plot(delay,tmp_d)
    
end

%% scatter of max locations
figure
subplot(4,1,[1:3])
hold on
col = summer(numel(MERGED));
params = [];
for iF=1:numel(MERGED)
    [~,max_loc]=max(MERGED(iF).average_triggered,[],2);
    tmp_d=bins;
    tmp_d=tmp_d-bins(MERGED(iF).max_depth);
    valid_idx = max_loc'>1 & tmp_d>-500 & tmp_d<500;

    plot(tvec_spikes(max_loc),tmp_d,'.','Color',col(iF,:))
    p=polyfit(tmp_d(valid_idx),tvec_spikes(max_loc(valid_idx)),1);
    params=cat(1,params,p);
    delay=polyval(p,tmp_d);
    %plot(delay,tmp_d)
end
av_delay = polyval(median(params),[-2000 3000]);
plot(av_delay,[-2000 3000],'k','LineWidth',2)
subplot(4,1,4)
speed = [-2000 1]*params'*1000/2;
boxplot(speed)
title(sprintf('Delay: %.2f (%.2f,%.2f)ms/mm',quantile(speed,[.5 .1 .9])))
%% find delay between max chan and rest
figure
hold on
for iF=1:numel(MERGED)
    template = MERGED(iF).average_triggered(MERGED(iF).max_depth,:);
    [dd,mc]=finddelay(template',MERGED(iF).average_triggered');
    tmp_d=bins;
    tmp_d=tmp_d-bins(MERGED(iF).max_depth);
    valid = mc>0.1;
    scatter3(tvec_spikes(dd(valid)+101),tmp_d(valid),mc(valid),'b.')
end
%% max loc and delay as a function of histology
figure
hold on
for iF= 1:numel(MERGED)
session = MERGED(iF).name;
session_parts = strsplit(session,'_');
animal = session_parts{1};
exp_date = session_parts{2};
%mod_idx=quantile(MERGED(iF).average_triggered,[.01 .99],2);
    %mod_idx=(mod_idx(:,2)-mod_idx(:,1))./mean(MERGED(iF).average_triggered,2);
    template = MERGED(iF).average_triggered(MERGED(iF).max_depth,:);
    [dd,mc]=finddelay(template',MERGED(iF).average_triggered');
table_range = strcmp(histology.animal,animal);
exp_row=all(histology.date==exp_date,2);
idx = find(table_range & exp_row);
if any(idx)
    
    pos3D = histology.origin(idx,:) + bins'*histology.unit_vector(idx,:);
    valid_idx = mc>0.1;
    cmap = winter(nnz(valid_idx));
    sz=linspace(1,15,nnz(valid_idx));
    [a,b]=sort(dd(valid_idx));
    [a,b_sz]=sort(abs(dd(valid_idx)));
    scatter3(pos3D(valid_idx,1),pos3D(valid_idx,2),pos3D(valid_idx,3),sz(b_sz),cmap(b,:))
end
end
xlabel('ML')
ylabel('Layer')
zlabel('Depth')
axis equal
%%
%% theta modulation idx, unalign, max chan, highest_chan, 'final depth'

zz=NaN(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
gg=NaN(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
hh=NaN(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
ii=hh;
n_bins = numel(bins);
d_vec = [-3840:40:3840];

theta_mat=zeros(97,numel(MERGED));
for iF=1:numel(MERGED)
    offset = 97-MERGED(iF).max_depth+1;
    template = MERGED(iF).average_triggered(MERGED(iF).max_depth,:);
    [dd,mc]=finddelay(template',MERGED(iF).average_triggered');
    mod_idx = dd;
    mod_idx(mc<0.7)=NaN;
    theta_mat(:,iF)=mod_idx;
    zz(offset:offset+n_bins-1,iF)=mod_idx;
    
    highest_chan = strfind(mean(MERGED(iF).average_triggered,2)'>0,[1 0 0 0 0])+1;
    if numel(highest_chan)>1
        if highest_chan(1)>=20
            highest_chan=highest_chan(1);
        else
            highest_chan=highest_chan(2);
        end
    end
    offset = 97-highest_chan+1;
    gg(offset:offset+n_bins-1,iF)=mod_idx;
    if ~isnan(final_depth(iF))
    offset = 97-lowest_depth(iF);
    hh(offset:offset+n_bins-1,iF)=mod_idx;
    end
    if ~isnan(mec_entry(iF))
    offset = 97-entry_bin(iF);
    ii(offset:offset+n_bins-1,iF)=mod_idx;
    end
end
%sortidx = 1:90;
[~,sortidx]=sort(ypos);
figure
subplot(1,5,1)
imagesc(flipud(theta_mat(:,sortidx)),[-10 10])
subplot(1,5,2)
imagesc(flipud(zz(:,sortidx)),[-10 10])
subplot(1,5,3)
imagesc(flipud(gg(:,sortidx)),[-10 10])
title('Brain Entry Act')
subplot(1,5,4)
imagesc(flipud(hh(:,sortidx)),[-10 10])
title('Brain Entry Histo')
subplot(1,5,5)
imagesc(flipud(ii(:,sortidx)),[-10 10])
title('MEC Entry')

figure
yl=col_range;
subplot(5,1,1)
plot(bins,nanmean(theta_mat,2))
ylim(yl)
subplot(5,1,2)
plot(d_vec,nanmean(zz,2))
ylim(yl)
subplot(5,1,3)
plot(d_vec,nanmean(gg,2))
ylim(yl)
subplot(5,1,4)
plot(d_vec+abs(min(d_vec)),nanmean(hh,2))
ylim(yl)
subplot(5,1,5)
plot(d_vec,nanmean(ii,2))
ylim(yl)



%%
figure
plot([-2000 3000],av_delay)
%% average spikemats, align to max chan

zz=zeros(2*size(MERGED(iF).average_triggered,1),size(MERGED(iF).average_triggered,2),length(MERGED));
n_bins = numel(bins);
for iF=1:length(MERGED)
    offset = 97-MERGED(iF).max_depth;
    aa_norm = bsxfun(@rdivide,MERGED(iF).average_triggered,sum(MERGED(iF).average_triggered,2));
    zz(offset:offset+n_bins-1,:,iF)=aa_norm;
end
figure
imagesc(flipud(nanmean(zz,3)),[0 0.007])

%set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
set(gca,'YTick',[1 50 99 99+50 194],'YTickLabel',[-max(bins) -max(bins)/2 0 max(bins)/2 max(bins)])
xlabel('Time [s]')
ylabel('Distance from max channel')
%%
histology = process_histology;
final_depth = nan(1,numel(MERGED));
mec_entry = final_depth;
ypos=final_depth;
xpos= final_depth;
for iF= 1:numel(MERGED)
session = MERGED(iF).name;
session_parts = strsplit(session,'_');
animal = session_parts{1};
exp_date = session_parts{2};

table_range = strcmp(histology.animal,animal);
exp_row=all(histology.date==exp_date,2);
idx = find(table_range & exp_row);
tmp_d=histology.FinalDepth(idx);
tmp_mec = histology.Z2(idx);
if ~isempty(tmp_d)
final_depth(iF)=tmp_d;
end
if ~isempty(tmp_mec) && tmp_mec>0
    mec_entry(iF)=tmp_mec;
    xpos(iF)=(histology.X1(idx)+histology.X2(idx))/2;
    ypos(iF)=histology.Y1(idx);
end
end
[aa,lowest_depth]=min(abs((bins-final_depth')),[],2);
[aa,entry_bin]=min(abs((bins-mec_entry')),[],2);
%%
figure
plot(bins(lowest_depth),final_depth,'.')

%% theta modulation idx, unalign, max chan, highest_chan, 'final depth'

mat_highest_channel=zeros(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
mat_activity=zeros(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
mat_insertion=zeros(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
mat_mec_entry=mat_insertion;
n_bins = numel(bins);
d_vec = [-3840:40:3840];

theta_mat=zeros(97,numel(MERGED));
for iF=1:numel(MERGED)
    offset = 97-MERGED(iF).max_depth+1;
    mod_idx=quantile(MERGED(iF).average_triggered,[.01 .99],2);
    mod_idx=(mod_idx(:,2)-mod_idx(:,1))./mean(MERGED(iF).average_triggered,2);
    theta_mat(:,iF)=mod_idx;
    mat_highest_channel(offset:offset+n_bins-1,iF)=mod_idx;
    
    highest_chan = strfind(mean(MERGED(iF).average_triggered,2)'>0,[1 0 0 0 0])+1;
    if numel(highest_chan)>1
        if highest_chan(1)>=20
            highest_chan=highest_chan(1);
        else
            highest_chan=highest_chan(2);
        end
    end
    offset = 97-highest_chan+1;
    mat_activity(offset:offset+n_bins-1,iF)=mod_idx;
    if ~isnan(final_depth(iF))
    offset = 97-lowest_depth(iF);
    mat_insertion(offset:offset+n_bins-1,iF)=mod_idx;
    end
    if ~isnan(mec_entry(iF))
    offset = 97-entry_bin(iF);
    mat_mec_entry(offset:offset+n_bins-1,iF)=mod_idx;
    end
end
%sortidx = 1:90;
col_range = [0 3];
[~,sortidx]=sort(ypos);
figure
subplot(1,5,1)
imagesc(flipud(theta_mat(:,sortidx)),col_range)
subplot(1,5,2)
imagesc(flipud(mat_highest_channel(:,sortidx)),col_range)
subplot(1,5,3)
imagesc(flipud(mat_activity(:,sortidx)),col_range)
title('Brain Entry Act')
subplot(1,5,4)
imagesc(flipud(mat_insertion(:,sortidx)),col_range)
title('Brain Entry Histo')
subplot(1,5,5)
imagesc(flipud(mat_mec_entry(:,sortidx)),col_range)
title('MEC Entry')

figure
yl=col_range;
subplot(5,1,1)
plot(bins,nanmean(theta_mat,2))
ylim(yl)
subplot(5,1,2)
plot(d_vec,nanmean(mat_highest_channel,2))
ylim(yl)
subplot(5,1,3)
plot(d_vec,nanmean(mat_activity,2))
ylim(yl)
subplot(5,1,4)
plot(d_vec+abs(min(d_vec)),nanmean(mat_insertion,2))
ylim(yl)
subplot(5,1,5)
plot(d_vec,nanmean(mat_mec_entry,2))
ylim(yl)



%% plot probe tracts
figure
hold on
for iF= 1:numel(MERGED)
session = MERGED(iF).name;
session_parts = strsplit(session,'_');
animal = session_parts{1};
exp_date = session_parts{2};
mod_idx=quantile(MERGED(iF).average_triggered,[.01 .99],2);
    mod_idx=(mod_idx(:,2)-mod_idx(:,1));%./mean(MERGED(iF).average_triggered,2);
table_range = strcmp(histology.animal,animal);
exp_row=all(histology.date==exp_date,2);
idx = find(table_range & exp_row);
if any(idx)
    
    pos3D = histology.origin(idx,:) + bins'*histology.unit_vector(idx,:);
    valid_idx = ~isnan(mod_idx);
    cmap = summer(nnz(valid_idx));
    sz=linspace(1,15,nnz(valid_idx));
    [a,b]=sort(mod_idx(valid_idx));
    scatter3(pos3D(valid_idx,1),pos3D(valid_idx,2),pos3D(valid_idx,3),sz(b),cmap(b,:))
end
end
xlabel('ML')
ylabel('Layer')
zlabel('Depth')
axis equal
%% scatter depth and ml
figure
hold on
for iF= 1:numel(MERGED)
session = MERGED(iF).name;
session_parts = strsplit(session,'_');
animal = session_parts{1};
exp_date = session_parts{2};
mod_idx=quantile(MERGED(iF).average_triggered,[.01 .99],2);
    mod_idx=(mod_idx(:,2)-mod_idx(:,1))./mean(MERGED(iF).average_triggered,2);
table_range = strcmp(histology.animal,animal);
exp_row=all(histology.date==exp_date,2);
idx = find(table_range & exp_row);
if any(idx)
    
    pos3D = histology.origin(idx,:) + bins'*histology.unit_vector(idx,:);
    valid_idx = ~isnan(mod_idx);
    cmap = summer(nnz(valid_idx));
    [a,b]=sort(mod_idx(valid_idx));
    subplot(1,2,1)
    hold on
    scatter(pos3D(valid_idx,2),pos3D(valid_idx,3),1,cmap(b,:))
        subplot(1,2,2)
    hold on
    scatter(pos3D(valid_idx,1),pos3D(valid_idx,3),1,cmap(b,:))
end
end
subplot(1,2,1)
xlabel('Layer')
ylabel('Depth')
axis equal
subplot(1,2,2)
xlabel('ML')
axis equal
ylabel('Depth')
