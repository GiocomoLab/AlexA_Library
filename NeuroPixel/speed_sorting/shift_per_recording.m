%% struct_names
dataset = 'shiftFilter_corrected';
path = ['Z:\giocomo\attialex\speed_' dataset];
filenames = dir(fullfile(path,'*.mat'));
ops = load(fullfile(path,'parameters'));
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
        if startsWith(r,'VISp')
            r='VISp';
        end
        try
            [ma,mi]=max(a.stability(iidx,:),[],2);
            stab = a.stability(iidx,26);
            tmp = [stab,mi,ones(size(ma))*iF];
            
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
        idx = data_this(:,3)==a(iP) & data_this(:,1)>.20;
        averages(iP)=nanmean(ops.factors(data_this(idx,2)));
    end
    AV_ALL.(reg{iR}) = averages;
    AV{iR}=averages;
end
%%
figure
data_this = STABILITY.VISp;
idx = data_this(:,1)>.2;
tmp = ops.factors(data_this(idx,2));
histogram(tmp,'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])

figure
plotSpread(AV)
set(gca,'XTickLabel',reg,'XTickLabelRotation',90)

%%
reg = {'VISp','RSC','MEC','RHP'};
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
AV={};
figure
for iR = 1:numel(reg)
    subplot(numel(reg),1,iR)
    data_this = STABILITY.(reg{iR});
    [a,b]=unique(data_this(:,3));
    averages  = zeros(size(a));
    for iP = 1:numel(a)
        idx = data_this(:,3)==a(iP) & data_this(:,1)>.20;
        averages(iP)=nanmean(ops.factors(data_this(idx,2)));
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
saveas(gcf,'F:/temp/figures/shift_raincloud.pdf')
%%
figure
plotSpread(AV)
set(gca,'XTickLabel',reg,'XTickLabelRotation',45)
ylabel('average speed correction lag')