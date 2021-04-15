%%
%data_table = readtable("C:\Users\giocomolab\Downloads\vi_trippy_pairs.xlsx");
data_table = readtable('/Volumes/T7/attialex/receptive_fields/vi_trippy_pairs.xlsx');
%trippy_path = 'F:\Alex\receptive_fields';
trippy_path = '/Volumes/T7/attialex/receptive_fields';
%projection = load('C:\code\AlexA_Library\NeuroPixel\receptive_fields\dorsal_projection.mat');
projection = load('/users/attialex/code/AlexA_Library/NeuroPixel/receptive_fields/dorsal_projection.mat');
projection = projection.projection;
uA= unique(projection);
annotation_volume_location = 'C:\code\allenCCF\Allen\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\code\allenCCF\Allen\structure_tree_safe_2017.csv';
template_volume_location = 'C:\code\allenCCF\Allen\template_volume_10um.npy';
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end
%%
groups=unique(data_table.group);
marker = {'x','*','o'};
cols = brewermap(3,'Set1');
for iG=1%:numel(groups)
    idx = find(data_table.group == groups(iG));
    animal = data_table.Animal(idx(1));
    probe_nbrs = data_table.Probe_nbr(idx);
    [~,pp]= plotProbeTrackAnimal(animal,'',av,st,probe_nbrs');
    if length(idx)<2
        continue
    end
    figure('Renderer','Painters','Color','white')
    for ii=1:numel(idx)
        data = load(fullfile(trippy_path,data_table.Trippy_name{idx(ii)},'receptive_fields.mat'));
        
        allMU=[];
        
        for iF=1:numel(data.fields)
            p=data.fields{iF};
            MU=[];
            for iROI=1:numel(p)
                if p{iROI}.field_sign ==-1
                    
                    col='b';
                else
                    col = 'r';
                end
%                 subplot(1,3,1)
%                 hold on
%                 plot(p{iROI}.xy(1:3:end,1),p{iROI}.xy(1:3:end,2),strcat(marker{ii},col))
                
                
                MU=cat(1,MU,p{iROI}.mu);
            end
            if size(MU,1)>1
                MU = mean(MU);
            end
            if ~isnan(MU)
                
                %scatter(MU(1),MU(2),15,cols(ii,:),marker{ii})
                allMU=cat(1,allMU,MU);
            end
        end
        subplot(1,2,1)
        hold on
        scatter(allMU(:,1),allMU(:,2),20,cols(ii,:),marker{ii})
        axis image
        set(gca,'YDir','reverse')

        subplot(1,2,2)
        for iA=2:numel(uA)
            BW = projection == uA(iA);
            [B,L] = bwboundaries(BW,'noholes');
            hold on
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1),'Color', [.2 .2 .2], 'LineWidth', 1)
            end
        end
        hold on
        hold on
        for aP=1:numel(pp)
            scatter(pp{aP}(2,3),pp{aP}(2,1),32,cols(aP,:),marker{aP})
        end
        
        bregma = allenCCFbregma();
        plot(bregma(3),bregma(1),'x')
        xline(bregma(3),'--')
        ylim([4        1318])
        xlim([ 54        1088])
        xlabel('L->R')
        ylabel('P->A')
        set(gca,'YDir','reverse')
        axis image
        xlim([600 990])
        ylim([600 1200])
        
    end
title('Probe entry points')

    subplot(1,2,1)
    title('RF Centers on screen')
    xlabel('lateral -> medial')
    ylabel('down -> up')
    xlim([0 176])
    ylim([0 96])
    %legend(data_table.Trippy_name(idx),'Interpreter','None')
end