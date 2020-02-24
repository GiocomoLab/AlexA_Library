%%
addpath(genpath('/home/users/attialex/AlexA_Library'));
addpath(genpath('/home/users/attialex/spikes/'));
%% find sessions based on table

%session_info=readtable('/oak/stanford/groups/giocomo/attialex/NP_DATA/session_summary.xlsx');

filenames = dir('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*dark*');

%%
root='/oak/stanford/groups/giocomo/attialex/NP_DATA';
im_save_dir= '/oak/stanford/groups/giocomo/attialex/wave_data_dark/images';
data_save_dir = '/oak/stanford/groups/giocomo/attialex/wave_data_dark/data';
if ~isfolder(im_save_dir)
    mkdir(im_save_dir);
end
if~isfolder(data_save_dir)
    mkdir(data_save_dir);
end
%root = 'Y:\giocomo\attialex\NP_DATA';
%root = fullfile('/oak/stanford/groups/giocomo','export','data','Projects','JohnKei_NPH3','E1','E1_190614_johncontrasttrack_train1_g0');
%Files = dir(fullfile(root,'*_contrast_*.mat'));
%Files = dir(fullfile(root,'*train*.mat'));

%%
p=gcp('nocreate');
if isempty(p)
    parpool(12);
end
%%
bins = 0:40:3840;

for iF=1:numel(filenames)

    fn = fullfile(root,filenames(iF).name);
    dataset = load(fn);
    [~,session_name] = fileparts(fn);
    try
        if isfile(fullfile('/oak/stanford/groups/giocomo/attialex/distance_coding3/data',[session_name '.mat']))
        acg=load(fullfile('/oak/stanford/groups/giocomo/attialex/distance_coding3/data',session_name));
        [aa_spikes,maxChan_spikes,spikefig,aa_mec_spikes]=quantification_session(dataset,bins,acg);
        else
        [aa_spikes,maxChan_spikes,spikefig,aa_mec_spikes]=quantification_session(dataset,bins);
        end
    set(spikefig,'Renderer','painters');
    saveas(spikefig,sprintf('%s/%s.pdf',im_save_dir,session_name))
    close(spikefig)
    matfi = matfile(fullfile(data_save_dir,session_name),'Writable',true);
    matfi.average_triggered = aa_spikes;
    matfi.mec_triggered = aa_mec_spikes;
    matfi.max_depth = maxChan_spikes;
    matfi.anatomy = dataset.anatomy;
    matfi.bins = bins;
    catch ME
        ME.message
        warning(strcat('something wrong with ',fn))
    end
    
    
end
%
%% scatter of max locations
time_bins = 0.002; 

tvec_spikes = [-100:100]*time_bins;

% figure
% hold on
% ANA_FILES = dir(fullfile(data_save_dir,'*.mat'));
% col = summer(numel(ANA_FILES));
% 
% for iF=1:numel(ANA_FILES)
%     
%     matfi = matfile(fullfile(ANA_FILES(iF).folder,ANA_FILES(iF).name));
%     [~,max_loc]=max(matfi.average_triggered,[],2);
%     tmp_d=bins;
%     tmp_d=tmp_d-bins(matfi.max_depth);
%     plot(tvec_spikes(max_loc),tmp_d,'.','Color',col(iF,:))
%         p=polyfit(tmp_d,tvec_spikes(max_loc),1);
%     delay=polyval(p,tmp_d);
%     plot(delay,tmp_d)
%     
% end

%% scatter of max locations
figure
subplot(4,2,[1 3 5])
hold on
ANA_FILES = dir(fullfile(data_save_dir,'*.mat'));
%ANA_FILES = dir(fullfile(im_save_dir,'good_waves','*.png'));

col = summer(numel(ANA_FILES));
params = [];
offsets=[];
for iF=1:numel(ANA_FILES)
        matfi = matfile(fullfile(data_save_dir,[ANA_FILES(iF).name(1:end-3) 'mat'] ));

    [~,max_loc]=max(matfi.average_triggered,[],2);
    %max_loc(max_loc==1)=[];
    bins = matfi.bins;
    tmp_d=bins;
    tmp_d=tmp_d-bins(matfi.max_depth);
    valid_idx = max_loc'>1 & tmp_d<400 & tmp_d>-1000;
    plot_idx = max_loc>1;
    plot(tvec_spikes(max_loc(plot_idx)),tmp_d(plot_idx),'.','Color',col(iF,:))
    p=polyfit(tmp_d(valid_idx),tvec_spikes(max_loc(valid_idx)),1);
    params=cat(1,params,p);
    delay=polyval(p,tmp_d);
    anatomy = matfi.anatomy;
    [~,mec_bin]=min(abs(bins-anatomy.z2));
    offsets(end+1)= anatomy.z2-bins(matfi.max_depth);
    %plot(delay,tmp_d)
end
av_delay = polyval(median(params),[-2000 2000]);
plot(av_delay,[-2000 2000],'k','LineWidth',2)

xlabel('Delay [ms]')
ylabel('Distance from peak channel')
yl = ylim();


subplot(4,2,7)
speed = [-2000 1]*params'*1000/2;
boxplot(speed)
title(sprintf('Delay: %.2f (%.2f,%.2f)ms/mm',quantile(speed,[.5 .1 .9])))
subplot(4,2,[2 4 6])
edges = -3000:500:3000;
counts = histcounts(offsets,edges);
cent = (edges(1:end-1)+edges(2:end))/2;
plot(counts,cent)
% subplot(4,1,[1:3])
% for iF=1:numel(offsets)
% plot(0,offsets(iF),'kx')
% end
ylim(yl);
plotSpread(offsets')
set(gcf,'Renderer','Painters')
saveas(gcf,'/oak/stanford/groups/giocomo/attialex/FIGURES/wave_1.pdf')
%% scatter of max locations relative to mec
figure
subplot(4,2,[1 3 5])
hold on
ANA_FILES = dir(fullfile(im_save_dir,'good_waves','*.png'));
col = summer(numel(ANA_FILES));
params = [];
offsets=[];
for iF=1:numel(ANA_FILES)
        matfi = matfile(fullfile(data_save_dir,[ANA_FILES(iF).name(1:end-3) 'mat']));

    [~,max_loc]=max(matfi.mec_triggered,[],2);
    %max_loc(max_loc==1)=[];
    bins = matfi.bins;
    tmp_d=bins;
  ana = matfi.anatomy;
    z2=ana.z2;
        [~,mec_bin]=min(abs(bins-z2));
        tmp_d = tmp_d-bins(mec_bin);
   
    
    if isnan(z2)
        continue
    end
    valid_idx = max_loc'>1 & tmp_d<-100;% & tmp_d>-1000;
    plot_idx = max_loc'>1;% & tmp_d<400; 
    plot(tvec_spikes(max_loc(plot_idx)),tmp_d(plot_idx),'.','Color',col(iF,:))
    p=polyfit(tmp_d(valid_idx),tvec_spikes(max_loc(valid_idx)),1);
    params=cat(1,params,p);
    delay=polyval(p,tmp_d);
    anatomy = matfi.anatomy;
    [~,mec_bin]=min(abs(bins-anatomy.z2));
    plot(delay,tmp_d)
end
av_delay = polyval(mean(params),[-2000 2000]);
plot(av_delay,[-2000 2000],'k','LineWidth',2)

xlabel('Delay [ms]')
ylabel('Distance from mec - 800')
yl = ylim();


subplot(4,2,7)
speed = [-2000 1]*params'*1000/2;
boxplot(speed)
title(sprintf('Delay: %.2f (%.2f,%.2f)ms/mm',quantile(speed,[.5 .1 .9])))
subplot(4,2,[2 4 6])

% subplot(4,1,[1:3])
% for iF=1:numel(offsets)
% plot(0,offsets(iF),'kx')
% end
ylim(yl);
set(gcf,'Renderer','Painters')
%saveas(gcf,'/oak/stanford/groups/giocomo/attialex/FIGURES/wave_max_mec.pdf')

%% make merged struct
ANA_FILES = dir(fullfile(data_save_dir,'*.mat'));
col = summer(numel(ANA_FILES));
MERGED = struct;
for iF=1:numel(ANA_FILES)
        matfi = matfile(fullfile(ANA_FILES(iF).folder,ANA_FILES(iF).name));
        fn = fieldnames(matfi);
        for ii=1:numel(fn)
            MERGED(iF).(fn{ii}) = matfi.(fn{ii});
        end
end
% 
% %% theta modulation idx, unalign, max chan, highest_chan, 'final depth'
% 
% zz=NaN(2*size(MERGED(iF).average_triggered,1)-1,numel(MERGED));
% 
% ii=zz;
% n_bins = numel(bins);
% d_vec = [-3840:40:3840];
% 
% theta_mat=zeros(97,numel(MERGED));
% for iF=1:numel(MERGED)
%     offset = 97-MERGED(iF).max_depth+1;
%     template = MERGED(iF).average_triggered(MERGED(iF).max_depth,:);
%     [dd,mc]=finddelay(template',MERGED(iF).average_triggered');
%     mod_idx = dd;
%     mod_idx(mc<0.7)=NaN;
%     theta_mat(:,iF)=mod_idx;
%     zz(offset:offset+n_bins-1,iF)=mod_idx;
%     
%     if ~isempty(MERGED(iF).anatomy) && ~isnan(MERGED(iF).anatomy.z2)
%         [aa,entry_bin]=min(abs((bins-MERGED(iF).anatomy.z2')));
%     offset = 97-entry_bin;
%     ii(offset:offset+n_bins-1,iF)=mod_idx;
%     end
% end
% %sortidx = 1:90;
% %[~,sortidx]=sort(ypos);
% sortidx = 1:size(theta_mat,2);
% figure
% subplot(1,3,1)
% imagesc(flipud(theta_mat(:,sortidx)),[-10 10])
% subplot(1,3,2)
% imagesc(flipud(zz(:,sortidx)),[-10 10])
% subplot(1,3,3)
% imagesc(flipud(ii(:,sortidx)),[-10 10])
% title('MEC Entry')
% 
% figure
% yl='auto'
% subplot(5,1,1)
% plot(bins,nanmean(theta_mat,2))
% ylim(yl)
% subplot(5,1,2)
% plot(d_vec,nanmean(zz,2))
% subplot(5,1,5)
% plot(d_vec,nanmean(ii,2))
% ylim(yl)




%% average spikemats, align to max chan

zz=zeros(2*size(MERGED(iF).average_triggered,1),size(MERGED(iF).average_triggered,2),length(MERGED));
n_bins = numel(bins);
offsets=[];
for iF=1:length(MERGED)
    offset = size(MERGED(iF).average_triggered,1)-MERGED(iF).max_depth;
    [~,mec_bin]=min(abs(bins-MERGED(iF).anatomy.z2));
    %offsets(end+1)= size(MERGED(iF).average_triggered,1)-mec_bin;
    offsets(end+1)=size(MERGED(iF).average_triggered,1)-(mec_bin-MERGED(iF).max_depth);
    aa_norm = bsxfun(@rdivide,MERGED(iF).average_triggered,sum(MERGED(iF).average_triggered,2));
    zz(offset:offset+n_bins-1,:,iF)=aa_norm;
end
figure
subplot(1,2,1)
tmp = flipud(nanmean(zz,3));
tmp = bsxfun(@rdivide,tmp,sum(tmp,2));
imagesc(tmp,[0 0.01])
colormap summer
%set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
set(gca,'YTick',[1 50 99 99+50 194],'YTickLabel',[-max(bins) -max(bins)/2 0 max(bins)/2 max(bins)])
xlabel('Time [s]')
ylabel('Distance from max channel')
hold on
%plot(100*ones(size(offsets)),offsets,'x')
subplot(1,2,2)
edges = 40:20:160;
counts = histcounts(offsets,edges);
cent = (edges(1:end-1)+edges(2:end))/2;
plot(counts,cent)
ylim([1 194])
set(gca,'YDir','reverse')
set(gcf,'Renderer','Painters')
saveas(gcf,'/oak/stanford/groups/giocomo/attialex/FIGURES/wave_2.pdf')
%
%%
zz=nan(2*size(MERGED(iF).average_triggered,1),size(MERGED(iF).average_triggered,2),length(MERGED));
n_bins = numel(bins);
for iF=1:length(MERGED)
    %offset = size(MERGED(iF).average_triggered,1)-MERGED(iF).max_depth;
    [~,mec_bin]=min(abs(bins-(MERGED(iF).anatomy.z2-800)));
    if isnan(MERGED(iF).anatomy.z2)
        continue
    end
    offset = size(MERGED(iF).average_triggered,1)-mec_bin;
    %aa_norm = bsxfun(@rdivide,MERGED(iF).mec_triggered,sum(MERGED(iF).mec_triggered,2));
    aa_norm = MERGED(iF).mec_triggered;
    zz(offset:offset+n_bins-1,:,iF)=zscore(aa_norm,0,2);
end
figure
imagesc(flipud(nanmean(zz,3)),[-1 1])

%set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
set(gca,'YTick',[1 50 99 99+50 194],'YTickLabel',[max(bins) max(bins)/2 0 -max(bins)/2 -max(bins)])
xlabel('Time [s]')
ylabel('Distance from MEC border -800')
set(gcf,'Renderer','Painters')
saveas(gcf,'/oak/stanford/groups/giocomo/attialex/FIGURES/wave_mec_zscore.pdf')
%
