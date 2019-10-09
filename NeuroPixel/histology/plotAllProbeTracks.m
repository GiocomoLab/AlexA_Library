% ------------------------------------------------------------------------
%          Display Probe Track
% ------------------------------------------------------------------------

%% ENTER PARAMETERS AND FILE LOCATION

% file location of probe points
%processed_images_folder = 'C:\Drive\Histology\for tutorial - sample data\Richards_done\processed';

% directory of reference atlas files

addpath(genpath('C:\code\npy-matlab'))
addpath(genpath('F:\code\allenCCF'))
% directory of histology

% name the saved probe points, to avoid overwriting another set of probes going in the same folder

% directory of reference atlas files
annotation_volume_location = 'F:\code\allenCCF\Allen\annotation_volume_10um_by_index.npy';
structure_tree_location = 'F:\code\allenCCF\Allen\structure_tree_safe_2017.csv';
template_volume_location = 'F:\code\allenCCF\Allen\template_volume_10um.npy';
afs = dir('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA*');
animals = {afs.name};





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

for iA = 1:numel(animals)
    plotProbeTrackAnimal(animals{iA},brainfig,av,st)
end