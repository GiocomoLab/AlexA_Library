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
data_table = readtable('C:\code\AlexA_Library\Neuropixel\histology\rsc_sessions.xlsx');
anatomy_table_loc = 'C:\code\AlexA_Library\Neuropixel\histology\rsc_tables';




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
for iS = 1:size(data_table,1)

    probe_number = data_table.Probe_nbr(iS);
    %[~,pp]= plotProbeTrackAnimal(data_table.Animal{iS},brainfig,av,st,probe_number);
    %AllProbePoints{end+1}=pp;
    [image_save_folder,probe_save_name_suffix,probe_lengths,processed_images_folder]=getProbeParametersAnimal(data_table.Animal{iS});
    border_table_path = fullfile(processed_images_folder,sprintf('probe_%d_region_table_new.mat',probe_number));
    borders_table = load(border_table_path);
    
    tbl_row.ProbeDepth = borders_table.probe_length/1000;
    data = load(fullfile('Z:\giocomo\attialex\NP_DATA_corrected',data_table.SessionName{iS}));
    metrics.tip_distance = data.anatomy.tip_distance;
    metrics.cluster_id = data.sp.cids;
    [cluster_table]=getClusterRegion_v2(tbl_row,borders_table.borders_table,borders_table.trajectory,metrics);
    %[cluster_region,cluster_parent,tip_distance,depth]=getClusterRegion(data.sp,probe_table.borders_table,probe_table.probe_length); % loock at _v2 to also give back x,y,z
    %turn into table
    
    %save table
    sn = [data_table.SessionName{iS} '_anatomy.csv'];
    fn = fullfile(anatomy_table_loc,sn);
    writetable(cluster_table,fn)
end
%%
