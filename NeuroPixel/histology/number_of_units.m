data_path = '/Volumes/Samsung_T5/attialex/NP_Data/';
files = dir(fullfile(data_path,'*.mat'));
animal_list = {};
valid_files = {};
for iF=1:numel(files)
    [~,sn]=fileparts(files(iF).name);
    
    parts = strsplit(sn,'_');
    animal = parts{1};
    date = parts{2};
    animal_date = [animal, '_',date];
    if ~ismember(animal_date,animal_list)
        animal_list{end+1}=animal_date;
        valid_files{end+1}=files(iF).name;
    end
end

%%
N_UNITS = struct();
SPAN = struct();
for iF=1:numel(valid_files)
    data = load(fullfile(data_path,valid_files{iF}));
    if ~isfield(data,'anatomy')
        sprintf('no anatomy for %s \n',valid_files{iF})
        continue
    end
    if isfield(data.anatomy,'parent_shifted')
        region = data.anatomy.parent_shifted;
    else
        region = data.anatomy.cluster_parent;
    end
    
    [a,b,c]=unique(region);
    good_idx = data.sp.cgs'==2;
    for iR = 1:numel(a)
        this_idx = c==iR & good_idx;
        n_units = nnz(this_idx);
        [l,s]=bounds(data.anatomy.tip_distance(this_idx));
        dd = abs(diff([l,s]));
        if nnz(this_idx)>1 && isvarname(a{iR})
            if isfield(N_UNITS,a{iR})
                N_UNITS.(a{iR}) = cat(1,N_UNITS.(a{iR}),n_units);
                SPAN.(a{iR}) = cat(1,SPAN.(a{iR}),dd);
            else
                N_UNITS.(a{iR})=n_units;
                SPAN.(a{iR}) = dd;
            end
        end
        
    end
end
%%

fn = fieldnames(N_UNITS);
N_UNITS.RSC = [];
SPAN.RSC = [];
for ifn = 1:numel(fn)
    if startsWith(fn{ifn},'RS')
        N_UNITS.RSC = cat(1,N_UNITS.RSC,N_UNITS.(fn{ifn}));
        SPAN.RSC = cat(1,SPAN.RSC,SPAN.(fn{ifn}));
    end
end
%%
fn = fieldnames(N_UNITS);
for ifn = 1:numel(fn)
    idx = N_UNITS.(fn{ifn})<10;
    N_UNITS.(fn{ifn})(idx)=[];
    SPAN.(fn{ifn})(idx)=[];
end




%%
% toplot = {'MEC','VISp','RSC'};
% TMP = struct();
% DENSITY = struct();
% X=[];
% G=[];
% for iV =1:numel(toplot)
%     TMP.(toplot{iV})=N_UNITS.(toplot{iV});
%     tmp = N_UNITS.(toplot{iV});%;./SPAN.(toplot{iV});
%     DENSITY.(toplot{iV}) = tmp;
%     X = cat(1,X,tmp);
%     G = cat(1,G,iV*ones(size(N_UNITS.(toplot{iV}))));
% end
% 
% figure('Position',[440   402   317   396])
% subplot(4,1,2)
% boxplot(X,G,'Symbol','k','Whisker',1,'Orientation','horizontal')
% lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% 
% % Change the boxplot color from blue to green
% a = get(get(gca,'children'),'children');   % Get the handles of all the objects
% %t = get(a,'tag');   % List the names of all the objects 
% %box1 = a(7);   % The 7th object is the first box
% set(a, 'Color', [0.2 0.2 0.2]);   % Set t
% set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
% b=plotSpread({DENSITY.MEC,DENSITY.VISp,DENSITY.RSC},'spreadWidth',.65,'xyOri','flipped','distributionMarkers','.','distributionColors','k');
% set(gca,'YTickLabel',regions)
% xlabel('# of Units')
% box off
% toplot = {'MEC','VISp','RSC'};
% TMP = struct();
% DENSITY = struct();
% X=[];
% G=[];
% for iV =1:numel(toplot)
%     TMP.(toplot{iV})=N_UNITS.(toplot{iV});
%     tmp = N_UNITS.(toplot{iV})./SPAN.(toplot{iV});
%     DENSITY.(toplot{iV}) = tmp;
%     X = cat(1,X,tmp);
%     G = cat(1,G,iV*ones(size(N_UNITS.(toplot{iV}))));
% end
% 
% subplot(4,1,1)
% boxplot(X,G,'Symbol','k','Whisker',1,'Orientation','horizontal')
% lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
% 
% % Change the boxplot color from blue to green
% a = get(get(gca,'children'),'children');   % Get the handles of all the objects
% %t = get(a,'tag');   % List the names of all the objects 
% %box1 = a(7);   % The 7th object is the first box
% set(a, 'Color', [0.2 0.2 0.2]);   % Set t
% b=plotSpread({DENSITY.MEC,DENSITY.VISp,DENSITY.RSC},'spreadWidth',.65,'xyOri','flipped','distributionMarkers','.','distributionColors','k');
% set(gca,'YTickLabel',regions)
% xlabel('Density [#units/um]')
% box off

%%
 cmap=cbrewer('qual','Set2',3,'pchip');
markers = {'o','x','v'};
toplot = {'MEC','VISp','RSC'};
TMP = struct();
DENSITY = struct();
X=[];
G=[];
figure('Position',[680   842   686   256])
for iV =1:numel(toplot)
    b_units=N_UNITS.(toplot{iV});
    span =SPAN.(toplot{iV});
    subplot(1,2,1)
    hold on
    scatter(b_units,span,25,cmap(iV,:),'o','filled')
    subplot(1,2,2)
    hold on
    scatter(b_units,span,25,[0 0 0],markers{iV})
end
subplot(1,2,1)
box off
legend(toplot)
legend('Location','southeast')
subplot(1,2,2)
box off
legend(toplot)
legend('Location','southeast')

saveas(gcf,'/Volumes/Samsung_T5/attialex/images/NUnitsvsSpan.pdf')


% %%
% figure
% subplot(2,2,1)
% violinplot(TMP)
% ylabel('Number of Units')
% box off
% subplot(2,2,2)
% h1 = raincloud_plot(TMP.VISp, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%     'box_col_match', 0,'cloud_edge_col', cb(1,:),'line_width',1);
% h2 = raincloud_plot(TMP.RSC, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'cloud_edge_col', cb(4,:),'line_width',1);
% h3 = raincloud_plot(TMP.MEC, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,...
%     'cloud_edge_col', cb(5,:),'line_width',1);
% legend('Visp','RSC','MEC');
% title('number of units');
% set(gca, 'YLim', [-8 15]*10e-4);
% 
% box off
% 
% 
% subplot(2,2,3)
% violinplot(DENSITY)
% ylabel('Units per um of probe')
% box off
% set(gcf,'Renderer','Painters')
% saveas(gcf,'/Volumes/Samsung_T5/attialex/images/NUnits.pdf')