%adapt table generating script to load anatomy struct from main data
sensory_delay = load('/Volumes/T7/attialex/sensory_delay.mat');
sensory_delay = sensory_delay.sensory_delay;
histo=process_histology();
%fill in depth, z2 and others from that
raw_data_path = '/Volumes/T7/attialex/NP_DATA_corrected';
%%
sessions = unique(sensory_delay.Session(startsWith(sensory_delay.BrainRegion,'MEC')));
extra_fields = {'z2','tip_distance','FinalDepth','x','y','z'};
if ~ismember(extra_fields{1},fieldnames(sensory_delay))
for iF=1:numel(extra_fields)
    sensory_delay.(extra_fields{iF}) = nan(size(sensory_delay.CellID));
end
for iS = 1:numel(sessions)
    if isfile(fullfile(raw_data_path,[sessions{iS} '.mat']))
    data=load(fullfile(raw_data_path,[sessions{iS} '.mat']),'anatomy','sp');
    else
        fprintf('%s not found \n',sessions{iS})
        continue
    end
    parts = strsplit(sessions{iS},'_');
    animal = parts{1};
    date = parts{2};
    idx_histo = strcmp(histo.date,date) & strcmp(histo.Mouse,animal);
    origin = histo.origin(idx_histo,:);
    unit_vec = histo.unit_vector(idx_histo,:);
    idx_delay = find(strcmp(sensory_delay.Session,sessions{iS}));
    for iC=idx_delay'
        cluID = sensory_delay.CellID(iC);
        sub_idx = data.sp.cids==cluID;
        sensory_delay.z2(iC)=data.anatomy.z2;
        sensory_delay.FinalDepth(iC)=data.anatomy.FinalDepth;
        sensory_delay.tip_distance(iC) = data.anatomy.tip_distance(sub_idx);
        if ~isnan(origin)
        pos = origin+data.anatomy.tip_distance(sub_idx)*unit_vec;
        sensory_delay.x(iC) = pos(1);
        sensory_delay.y(iC) = pos(2);
        sensory_delay.z(iC) = pos(3);
        end
    end
        
end
end
%%
sessions = unique(sensory_delay.Session(startsWith(sensory_delay.BrainRegion,'MEC')));
figure
hold on
for iS = 1:numel(sessions)
    
        idx_delay = find(strcmp(sensory_delay.Session,sessions{iS}));
        C=sensory_delay.DelayFactor(idx_delay)*-1;
    scatter3(sensory_delay.x(idx_delay),sensory_delay.y(idx_delay),sensory_delay.z(idx_delay),15,C,'filled')
end

xlabel('ML')
ylabel('Layer')
zlabel('Depth')
set(gca,'CLim',[-.1 .2])
axis equal
%%
aid = unique(sensory_delay.Session(startsWith(sensory_delay.BrainRegion,'MEC')));
midpoints= -2300:200:1500;
half_width = 200;
z2_temp = sensory_delay.z2;
%z2_temp(isnan(z2_temp))=1600;
DEPTH = sensory_delay.tip_distance-sensory_delay.z2;

allSites = nan(numel(aid),numel(midpoints));
for ii=1:numel(aid)
    valid_idx = strcmp(sensory_delay.Session,aid{ii});
    vals = nan(1,numel(midpoints));
    for iv=1:numel(vals)
        idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
        vals(iv)=nanmean(sensory_delay.DelayFactor(idx)*-1);
    end
    %plot(vals,edges(1:end-1))
    allSites(ii,:)=vals;
end




figure
subplot(1,2,1)
%plot(nanmean(allSites),edges(1:end-1))
err=nan(size(allSites,2),1);
for iS=1:numel(err)
    err(iS)=nanstd(allSites(:,iS));
    err(iS)=err(iS)/sqrt(nnz(~isnan(allSites(:,iS))));
end
%err = nanstd(allSites)/sqrt(numel(aid));
errorbar(nanmean(allSites),midpoints,[],[],err,err)
ylabel('Depth')
xlabel('Delay Factor')
grid on




DEPTH = sensory_delay.tip_distance-sensory_delay.z2;
valid_idx = startsWith(sensory_delay.Session,'np') ;
vals = nan(1,numel(midpoints));
for iv=1:numel(vals)
    idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    vals(iv)=nanmean(sensory_delay.DelayFactor(idx)*-1);
    err(iv)=nanstd(sensory_delay.DelayFactor(idx));
end
%plot(vals,edges(1:end-1))
hold on
%plot(vals,midpoints)
subplot(1,2,2)
hold on
for ii=1:size(allSites,1)
    plot(allSites(ii,:),midpoints)
end

%%
%%
dist = load('/Volumes/T7/attialex/dist_tuning_autocorr_bin2_smooth4.mat');
dist_tuning = dist.dist_tuning;
%%
d_tuned = dist_tuning.prom>=0.1 & dist_tuning.pval<=0.01;
dist_tuned = false(size(sensory_delay,1),1);
dark_data = dist_tuned;
for iC=1:size(sensory_delay,1)
    mouse = sensory_delay.Mouse{iC};
    date = sensory_delay.Date{iC};
    cellID = sensory_delay.CellID(iC);
    mouse_date = strcat(mouse,'_',date);
    
    idx = strcmp(dist_tuning.MouseDate,mouse_date) & dist_tuning.CellID == cellID;
    if nnz(idx)==1
    dist_tuned(iC)=d_tuned(idx);
    dark_data(iC)=true;
    elseif nnz(idx)>1
        idx = find(idx);
        dist_tuned(iC)=d_tuned(idx(1));
        dark_data(iC)=true;
    end
end
%%
%% plot delay factors for sessions with dark data only
aid = unique(sensory_delay.Session(startsWith(sensory_delay.BrainRegion,'MEC')));
midpoints= -2100:200:1500;
half_width = 200;
z2_temp = sensory_delay.z2;
%z2_temp(isnan(z2_temp))=1600;
DEPTH = sensory_delay.tip_distance-sensory_delay.z2;

allSitesTuned = nan(numel(aid),numel(midpoints));
allSitesNotTuned = allSitesTuned;
for ii=1:numel(aid)
    valid_idx = strcmp(sensory_delay.Session,aid{ii}) & dark_data;
    vals_tuned = nan(1,numel(midpoints));
    vals_notTuned = vals_tuned;
    for iv=1:numel(midpoints)
        idx_depth = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
        idx_tuned = idx_depth & dist_tuned;
        idx_nottuned = idx_depth & ~dist_tuned;
        vals_tuned(iv) = nanmean(sensory_delay.DelayFactor(idx_tuned)*-1);
        vals_notTuned(iv) = nanmean(sensory_delay.DelayFactor(idx_nottuned)*-1);
    end
    %plot(vals,edges(1:end-1))
    allSitesTuned(ii,:)=vals_tuned;
    allSitesNotTuned(ii,:) = vals_notTuned;
    
end

p=sum(~isnan(allSitesTuned),1);
allSitesTuned(:,p<2)=nan;

f = sum(~isnan(allSitesNotTuned),1);
allSitesNotTuned(:,f<2)=nan;

figure('Position',[440   133   410   665])
%plot(nanmean(allSites),edges(1:end-1))
hold on

errTuned=nan(size(allSitesTuned,2),1);
errNotTuned = errTuned;
diffTuned = allSitesTuned-allSitesNotTuned;
p_vals = errTuned;
for iS=1:numel(errTuned)
    errTuned(iS)=nanstd(allSitesTuned(:,iS));
    errTuned(iS)=errTuned(iS)/sqrt(nnz(~isnan(allSitesTuned(:,iS))));
    errNotTuned(iS)=nanstd(allSitesNotTuned(:,iS));
    errNotTuned(iS)=errNotTuned(iS)/sqrt(nnz(~isnan(allSitesNotTuned(:,iS))));
    try
    p_vals(iS)=signrank(diffTuned(:,iS));
    catch ME
    end
    
end
%err = nanstd(allSites)/sqrt(numel(aid));
errorbar(nanmean(allSitesTuned),midpoints,[],[],errTuned,errTuned)
hold on
errorbar(nanmean(allSitesNotTuned),midpoints,[],[],errNotTuned,errNotTuned)

for ip=1:numel(p_vals)
    if ~isnan(p_vals(ip));
        text(-.3,midpoints(ip),sprintf('%.3e',p_vals(ip)))
        
    end
end


ylabel('Depth')
xlabel('Delay Factor')
xlim([-.3 .3])
grid on
legend({'Tuned','Not Tuned'})
saveas(gcf,'/Users/attialex/Dropbox/temporary_images/delay_dark_tuned_depth.pdf')




%% plot delay factors for sessions with dark data only
aid = unique(sensory_delay.Session(startsWith(sensory_delay.BrainRegion,'MEC')));
midpoints= -2300:200:1500;
half_width = 200;
z2_temp = sensory_delay.z2;
%z2_temp(isnan(z2_temp))=1600;
DEPTH = sensory_delay.tip_distance-sensory_delay.z2;

allSites = nan(numel(aid),numel(midpoints));
for ii=1:numel(aid)
    valid_idx = strcmp(sensory_delay.Session,aid{ii}) & dark_data;
    vals = nan(1,numel(midpoints));
    for iv=1:numel(vals)
        idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
        vals(iv)=nanmean(sensory_delay.DelayFactor(idx)*-1);
    end
    %plot(vals,edges(1:end-1))
    allSites(ii,:)=vals;
end




figure
subplot(1,2,1)
%plot(nanmean(allSites),edges(1:end-1))
hold on

err=nan(size(allSites,2),1);
for iS=1:numel(err)
    err(iS)=nanstd(allSites(:,iS));
    err(iS)=err(iS)/sqrt(nnz(~isnan(allSites(:,iS))));
end
%err = nanstd(allSites)/sqrt(numel(aid));
errorbar(nanmean(allSites),midpoints,[],[],err,err)
ylabel('Depth')
xlabel('Delay Factor')
grid on



DEPTH = sensory_delay.tip_distance-sensory_delay.z2;
valid_idx = startsWith(sensory_delay.Session,'np') & dist_tuned;
vals = nan(1,numel(midpoints));
for iv=1:numel(vals)
    idx = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    vals(iv)=nanmean(sensory_delay.DelayFactor(idx)*-1);
    err(iv)=nanstd(sensory_delay.DelayFactor(idx));
end
plot(vals,midpoints)


legend({'all cells, averaged across sites','dist tuned cells only'})

subplot(1,2,2)
hold on

DEPTH = sensory_delay.tip_distance-sensory_delay.z2;
valid_idx = startsWith(sensory_delay.Session,'np') & dark_data;
vals_dt = nan(1,numel(midpoints));
vals_ndt = vals_dt;
err_tuned = vals_dt;
err_ndt = vals_dt;
p_vals = vals_ndt;
for iv=1:numel(vals_dt)
    idx_depth = DEPTH>midpoints(iv)-half_width & DEPTH<midpoints(iv)+half_width & valid_idx;
    idx_tuned = idx_depth & dist_tuned;
    idx_nt = idx_depth & ~dist_tuned;
    vals_dt(iv)=nanmedian(sensory_delay.DelayFactor(idx_tuned)*-1);
    err_tuned(iv)=nanstd(sensory_delay.DelayFactor(idx_tuned))/sqrt(nnz(~isnan(sensory_delay.DelayFactor(idx_tuned))));
    
    vals_ndt(iv)=nanmedian(sensory_delay.DelayFactor(idx_nt)*-1);
    err_ndt(iv)=nanstd(sensory_delay.DelayFactor(idx_nt))/sqrt(nnz(~isnan(sensory_delay.DelayFactor(idx_nt))));
    [~,p]=ttest2(sensory_delay.DelayFactor(idx_nt),sensory_delay.DelayFactor(idx_tuned));
    p_vals(iv)=p;
end
%plot(vals_dt,midpoints)
errorbar(vals_dt,midpoints,[],[],err_tuned,err_tuned)
errorbar(vals_ndt,midpoints,[],[],err_ndt,err_ndt)

for ip=1:numel(p_vals)
    if ~isnan(p_vals(ip));
        text(-.2,midpoints(ip),sprintf('%.3e',p_vals(ip)))
        
    end
end



grid on
legend({'dist tuned','not dist tuned'})

