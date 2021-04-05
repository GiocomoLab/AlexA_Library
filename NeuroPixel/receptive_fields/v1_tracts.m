% ------------------------------------------------------------------------
%          Display Probe Track
% ------------------------------------------------------------------------

%% ENTER PARAMETERS AND FILE LOCATION

% file location of probe points
%processed_images_folder = 'C:\Drive\Histology\for tutorial - sample data\Richards_done\processed';

% directory of reference atlas files

addpath(genpath('C:\code\npy-matlab'))
addpath(genpath('C:\code\allenCCF'))
% directory of histology

% name the saved probe points, to avoid overwriting another set of probes going in the same folder

% directory of reference atlas files
annotation_volume_location = 'C:\code\allenCCF\Allen\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\code\allenCCF\Allen\structure_tree_safe_2017.csv';
template_volume_location = 'C:\code\allenCCF\Allen\template_volume_10um.npy';
%afs = dir('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA*');
%animals = {afs.name};
data_table = readtable('G:\attialex\receptive_fields\vi_trippy_pairs.xlsx');
animals = unique(data_table.Animal);
animals = animals(2:end) %because first is empty




%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end



black_brain = true;
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';
brainfig = gcf;
AllProbePoints = {};
for iA = 1%:numel(animals)
   idx = startsWith(data_table.Animal,animals{iA});
   probe_numbers = data_table.Probe_nbr(idx);
   [~,pp]= plotProbeTrackAnimal(animals{iA},brainfig,av,st,probe_numbers');
   AllProbePoints{end+1}=pp;
end
%%
figure
for ii=1:5
    hold on
    plot3(pp{ii}(:,1),pp{ii}(:,2),pp{ii}(:,3))
end
xlim([4        1318])
zlim([ 54        1088])
ylim([17   755])
%%
figure
hold on
for iP=1:numel(AllProbePoints)
pp=AllProbePoints{iP};
for ii=1:numel(pp)
scatter(pp{ii}(2,3),pp{ii}(2,1))
end
end
bregma = allenCCFbregma()
plot(bregma(3),bregma(1),'x')
xline(bregma(3),'--')
ylim([4        1318])
xlim([ 54        1088])
xlabel('L->R')
ylabel('P->A')
set(gca,'YDir','reverse')
