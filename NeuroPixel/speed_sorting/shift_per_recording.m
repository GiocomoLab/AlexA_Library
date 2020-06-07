%% struct_names
dataset = 'filtered_new_11binfilt';
%dataset = 'filtered_new_5binfiltSpace_1binfiltSpeed';
%path = ['Z:\giocomo\attialex\speed_' dataset];
%path = '/Volumes/Samsung_T5/speed_filtered_greedy3';
path = '/Volumes/Samsung_T5/speed_filtered_new_22binspace_5binspeed2';
%path = '/Volumes/Samsung_T5/speed_filtered_greedy4';
fig_path ='/Volumes/Samsung_T5/attialex/images';
filenames = dir(fullfile(path,'*.mat'));
ops = load(fullfile(path,'parameters.mat'));
ops = ops.ops;
%% 

STABILITY = struct();

for iF=1:numel(filenames)
    try
    filename=filenames(iF).name;
    a=load(fullfile(path,filename));
    [unique_regions,~,ridx]=unique(a.region);
    for iR=1:numel(unique_regions)
        
        r=unique_regions{iR};
        iidx = ridx==iR;
        if startsWith(r,'RS')
            r='RSC';
        end
        if startsWith(r,'VISpm')
            r='VISm';
        end
        try
            dat = nanmean(a.all_factors(:,iidx))';
            
            stab = nanmean(a.all_stability(:,iidx))';
            fr_tmp = nanmean(a.FR(:,iidx))';
            tmp = [stab,dat,ones(size(dat))*iF,fr_tmp]; %stability, factors,recording_id, firing rate
            
        if ismember(r,fieldnames(STABILITY))
        STABILITY.(r) = cat(1,STABILITY.(r),tmp);
        else
            STABILITY.(r)=tmp;
        end
        catch ME
            disp(ME.message)
        end
    end
    catch ME
        disp(ME.message)
    end
        
end

%%
reg = fieldnames(STABILITY);
AV={};
for iR = 1:numel(reg)
    data_this = STABILITY.(reg{iR});
    [a,b]=unique(data_this(:,3));
    averages  = zeros(size(a));
    for iP = 1:numel(a)
        idx = data_this(:,3)==a(iP) & data_this(:,1)>.20 & data_this(:,4)>2;
        averages(iP)=nanmean((data_this(idx,2)));
    end
    AV_ALL.(reg{iR}) = averages;
    AV{iR}=averages;
end
%%
figure('Position',[1         735        1647         363],'Renderer','Painters')
regions = {'MEC','VISp','RSC'};
X=[];
G=[];
USED = struct();
TE = '';
TE2='';
for iR=1:numel(regions)
data_this = STABILITY.(regions{iR});
idx = data_this(:,1)>.1 & data_this(:,4)>1;
tmp = (data_this(idx,2));
subplot(1,3,1)
hold on
histogram(tmp,'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])
subplot(1,3,2)
hold on
ksdensity(tmp)
USED.(regions{iR})=tmp;
X=cat(1,X,tmp);
G=cat(1,G,iR*ones(size(tmp)));
TE=strcat(TE,sprintf('%s: %d of %d, ',regions{iR},nnz(idx),numel(idx)));
TE2=strcat(TE2,sprintf('%s: %.2f, +- %.3f. ',regions{iR},mean(tmp),std(tmp)/sqrt(nnz(idx))));
end
subplot(1,3,1)
text(-.3,.1,TE)
text(-.3,.11,TE2)
subplot(1,3,3)
violinplot(USED)
p1=ranksum(X(G==1),X(G==2));
p2 = ranksum(X(G==1),X(G==3));
text(2,0.2,sprintf('%3e',p1))
text(3,0.25,sprintf('%3e',p2))
saveas(gcf,'/Volumes/Samsung_T5/attialex/shift_all_cells.pdf')
% figure
% plotSpread(AV)
% set(gca,'XTickLabel',reg,'XTickLabelRotation',90)

%%
reg = {'VISp','RSC','MEC','RHP','ECT'};
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
AV={};
AV_ALL={};
figure
for iR = 1:numel(reg)
    subplot(numel(reg),1,iR)
    data_this = STABILITY.(reg{iR});
    %data_this(data_this<=-.15)=nan;
    [a,b]=unique(data_this(:,3));
    averages  = zeros(size(a));
    for iP = 1:numel(a)
        idx = data_this(:,3)==a(iP) & data_this(:,1)>.10 & data_this(:,4)>1;
        averages(iP)=nanmean((data_this(idx,2)));
        if nnz(idx)<nnz(data_this(:,3)==a(iP))*.1
            averages(iP)=nan;
        end
    end
    AV_ALL.(reg{iR}) = averages;
    AV{iR}=averages;
    raincloud_plot(averages, 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount',...
.3, 'dot_dodge_amount', .3, 'color', cb(1,:), 'cloud_edge_col', cb(1,:));
title(reg{iR});
set(gca, 'XLim', [-.2 .2]);
box off
end

figure
h1 = raincloud_plot(AV{1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'cloud_edge_col', cb(1,:),'line_width',1);
h2 = raincloud_plot(AV{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'cloud_edge_col', cb(4,:),'line_width',1);
 h3 = raincloud_plot(AV{3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,...
     'cloud_edge_col', cb(5,:),'line_width',1);
legend([h1{1} h2{1} h3{1}], reg);
title('Average shift per recording');
set(gca,'XLim', [-.25 .15], 'YLim', [-20 31]);

box off
set(gcf,'Renderer','Painters')
%saveas(gcf,'F:/temp/figures/shift_raincloud.pdf')

figure
plotSpread(AV)
set(gca,'XTickLabel',reg,'XTickLabelRotation',45)
ylabel('average speed correction lag')
%%

% %% Firing rate
% reg = {'VISp','RSC','MEC','RHP','ECT'};
% [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
% AV={};
% AV_ALL={};
% 
% for iR = 1:numel(reg)
%     data_this = STABILITY.(reg{iR});
%     [a,b]=unique(data_this(:,3));
%     averages  = zeros(size(a));
%     for iP = 1:numel(a)
%         idx = data_this(:,3)==a(iP) & data_this(:,1)>.1;
%         averages(iP)=nanmean((data_this(idx,4)));
%         if nnz(idx)<nnz(data_this(:,3)==a(iP))*.1
%             averages(iP)=nan;
%         end
%     end
%     AV_ALL.(reg{iR}) = averages;
%     AV{iR}=averages;
%     
% end
% 
% figure
% h1 = raincloud_plot(AV{1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
%      'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%      'box_col_match', 0,'cloud_edge_col', cb(1,:),'line_width',1);
% h2 = raincloud_plot(AV{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
%      'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'cloud_edge_col', cb(4,:),'line_width',1);
%  h3 = raincloud_plot(AV{3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
%      'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,...
%      'cloud_edge_col', cb(5,:),'line_width',1);
% legend([h1{1} h2{1} h3{1}], reg);
% title('Average Firing Rate');
% %set(gca,'XLim', [-.25 .15], 'YLim', [-20 31]);
% 
% box off
% set(gcf,'Renderer','Painters')
% saveas(gcf,fullfile(fig_path,'/firingrate_raincloud.pdf'))
% 
% figure
% plotSpread(AV)
% set(gca,'XTickLabel',reg,'XTickLabelRotation',45)
% ylabel('firing rate')
% %%
% %% Stability
% reg = {'VISp','RSC','MEC','RHP','ECT'};
% [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
% AV={};
% AV_ALL={};
% 
% for iR = 1:numel(reg)
%     data_this = STABILITY.(reg{iR});
%     [a,b]=unique(data_this(:,3));
%     averages  = zeros(size(a));
%     for iP = 1:numel(a)
%         idx = data_this(:,3)==a(iP) & data_this(:,1)>.1;
%         averages(iP)=nanmean((data_this(idx,1)));
%         if nnz(idx)<nnz(data_this(:,3)==a(iP))*.1
%             averages(iP)=nan;
%         end
%     end
%     AV_ALL.(reg{iR}) = averages;
%     AV{iR}=averages;
%     
% end
% 
% figure
% h1 = raincloud_plot(AV{1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
%      'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%      'box_col_match', 0,'cloud_edge_col', cb(1,:),'line_width',1);
% h2 = raincloud_plot(AV{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
%      'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'cloud_edge_col', cb(4,:),'line_width',1);
%  h3 = raincloud_plot(AV{3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
%      'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,...
%      'cloud_edge_col', cb(5,:),'line_width',1);
% legend([h1{1} h2{1} h3{1}], reg);
% title('Average Spatial Correlation');
% set(gca,'XLim', [-.10 .7], 'YLim', [-4 11]);
% 
% box off
% set(gcf,'Renderer','Painters')
% saveas(gcf,fullfile(fig_path,'/stability_raincloud.pdf'))
% 
% figure
% plotSpread(AV)
% set(gca,'XTickLabel',reg,'XTickLabelRotation',45)
% ylabel('firing rate')