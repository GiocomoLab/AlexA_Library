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
    cntr=cntr+1;
%     if mod(iF,7) ~=0
%         close(spikefig)
%     end
    catch ME
        ME.message
        warning(strcat('something wrong with ',Files(iF).name))
    end

end
save(fullfile(root,'MERGED_FIRSTSESSIONS'),'MERGED');
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


%%
figure
plot([-2000 3000],av_delay)
%% average spikemats

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