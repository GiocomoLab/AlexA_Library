%% struct_names
dataset = 'filtered_new_11binfilt';
path = ['Z:\giocomo\attialex\speed_' dataset];
filenames = dir(fullfile(path,'*.mat'));

%% 

FACTOR = struct();
DEPTH = struct();
for iF=1:numel(filenames)
        filename=filenames(iF).name;

    if startsWith(filename,'para')
        continue
    end
    a=load(fullfile(path,filename));
    [unique_regions,~,ridx]=unique(a.region);
    for iR=1:numel(unique_regions)
        try
        r=unique_regions{iR};
        iidx = ridx==iR;
        factors = nanmean(a.all_factors)';
        if startsWith(r,'RS')
            r='RSC';
        end
        if startsWith(r,'VISp')
            r='VISp';
        end
        try
        if ismember(r,fieldnames(STABILITY))
        STABILITY.(r) = cat(1,STABILITY.(r),factors(iidx));
        else
            STABILITY.(r)=factors(iidx);
        end
        catch ME
            disp(ME.message)
        end
    
    catch ME
        disp(ME.message)
        end
    end
        
end
            
%%
ops = load(fullfile(path,'parameters'));
ops = ops.ops;
fn = fieldnames(STABILITY);
fn = {'RSC'};
histfig=figure;
for iF=1:numel(fn)
    dat = STABILITY.(fn{iF});
    
    figure(histfig)
    subplot(1,2,1)
    hold on
    histogram(dat,'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])
    subplot(1,2,2)
    hold on
    [f,xi]=ksdensity(dat);
    plot(xi,f)
end
    legend(fn)
    xlabel('speed correction lag [s]')
    
%%
fn = fieldnames(FACTOR);

N = zeros(numel(fn),1);
for iF=1:numel(fn)
    N(iF)=size(FACTOR.(fn{iF}),2);
end

[a,sid]=sort(N,'descend');
figure
for iF=1:5
    
    dat = FACTOR.(fn{sid(iF)});
    subplot(5,2,iF*2)
    idx = dat(2,:)>.1;
   
%     %idx = ma>0.01;
%      idx = dat(:,26)>0.2;
%     plot(ops.factors(mi),ma,'.')
    subplot(5,2,iF*2-1)
    hold on
    histogram(dat(1,idx),'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])
        legend(fn{sid(iF)})
    xlabel('speed correction lag [s]')
end
%%
[cb] = cbrewer('qual', 'Set3', 5, 'pchip');

fn = fieldnames(FACTOR);

N = zeros(numel(fn),1);
for iF=1:numel(fn)
    N(iF)=size(FACTOR.(fn{iF}),2);
end

[a,sid]=sort(N,'descend');
figure('Position',[680   188   290   790])
for iF=1:5
    
    dat = FACTOR.(fn{sid(iF)});
    
    subplot(5,1,iF)
    idx = dat(2,:)>.2;

    raincloud_plot(dat(1,idx), 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount',...
.3, 'dot_dodge_amount', .3, 'color', cb(iF,:), 'cloud_edge_col', cb(iF,:),'line_width',1);
title(fn{sid(iF)});
set(gca, 'XLim', [-.2 .2]);
box off
    
end
set(gcf,'Renderer','Painters')
%saveas(gcf,'F:/temp/figures/shift_raincloud_topfive.pdf')
    %%
    ops = load(fullfile(path,'parameters'));
ops = ops.ops;
fn = fieldnames(FACTOR);
fn = {'MEC','VISp','RSC'};
histfig=figure;
DAT={};
for iF=1:numel(fn)
    dat = FACTOR.(fn{iF});
    
    idx = dat(2,:)>0.1;
    DAT{iF}=dat(1,idx);
end
  
    figure
h1 = raincloud_plot(DAT{1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'cloud_edge_col', cb(1,:),'line_width',1);
h2 = raincloud_plot(DAT{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'cloud_edge_col', cb(4,:),'line_width',1);
 h3 = raincloud_plot(DAT{3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,...
     'cloud_edge_col', cb(5,:),'line_width',1);
legend([h1{1} h2{1} h3{1}], fn);
title('Average shift per recording');
set(gca,'XLim', [-.25 .15], 'YLim', [-20 31]);
set(gcf,'Renderer','Painters')
%saveas(gcf,'F:/temp/figures/shift_raincloud_allcells.pdf')
