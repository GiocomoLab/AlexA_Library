function [brainfig,ProbePoints] = plotProbeTrackAnimal(animalName,brainfig,annotated_volume,structure_tree,probe_numbers)

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

% name of the saved probe points
%probe_save_name_suffix = 'electrode_track_1';
[image_save_folder,probe_save_name_suffix,probe_lengths,processed_images_folder]=getProbeParametersAnimal(animalName);

% either set to 'all' or a list of indices from the clicked probes in this file, e.g. [2,3]
probes_to_analyze = 'all';  % [1 2]
if isempty(probe_numbers)
    probes_to_analyze = 'all';  % [1 2]
else
    probes_to_analyze = probe_numbers;
end

% -----------
% parameters
% -----------
% how far into the brain did you go from the surface, either for each probe or just one number for all -- in mm

% from the bottom tip, how much of the probe contained recording sites -- in mm
active_probe_length = 3.84;

% distance queried for confidence metric -- in um
probe_radius = 100;

% overlay the distance between parent regions in gray (this takes a while)
show_parent_category = false;

% plot this far or to the bottom of the brain, whichever is shorter -- in mm
distance_past_tip_to_plot = .5;

% set scaling e.g. based on lining up the ephys with the atlas
% set to *false* to get scaling automatically from the clicked points
scaling_factor = false;

% show a table of regions that the probe goes through, in the console
show_region_table = true;

% black brain?
black_brain = true;


% close all




%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE
if ~isempty(annotated_volume)
    av=annotated_volume;
end
if ~isempty(structure_tree)
    st=structure_tree;
end
% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

% load probe points
probePoints = load(fullfile(processed_images_folder, ['probe_points' probe_save_name_suffix]));
ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0];
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};
fwireframe = [];

% scale active_probe_length appropriately
active_probe_length = active_probe_length*100;

% determine which probes to analyze
if strcmp(probes_to_analyze,'all')
    probes = 1:size(probePoints.pointList.pointList,1);
    if numel(probes)==1
        probes = 1:size(probePoints.pointList.pointList,2);
    end
else
    probes = probes_to_analyze;
end





%% PLOT EACH PROBE -- FIRST FIND ITS TRAJECTORY IN REFERENCE SPACE
ProbePoints = {};
for selected_probe = probes
    
    % get the probe points for the currently analyzed probe
    try
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [3 2 1]);
    catch
        continue
    end
    
    % get user-defined probe length from experiment
    if length(probe_lengths) > 1
        probe_length = probe_lengths(selected_probe);
    else
        probe_length = probe_lengths;
    end
    
    % get the scaling-factor method to use
    if scaling_factor
        use_tip_to_get_reference_probe_length = false;
        reference_probe_length = probe_length * scaling_factor;
        disp(['probe scaling of ' num2str(scaling_factor) ' determined by user input']);
    else
        use_tip_to_get_reference_probe_length = true;
        disp(['getting probe scaling from histology data...']);
    end
    
    % get line of best fit through points
    % m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
    [start_point,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
    if isnan(start_point(1))
        disp(['no points found for probe ' num2str(selected_probe)])
        continue
    end
    
    % ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
    if p(2)<0
        p = -p;
    end
    
    % determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
    ann = 10;
    isoCtxId = num2str(st.id(strcmp(st.acronym, 'Isocortex')));
    gotToCtx = false;
    while ~(ann==1 && gotToCtx)
        start_point = start_point-p; % step 10um, backwards up the track
        ann = av(round(start_point(1)),round(start_point(2)),round(start_point(3))); %until hitting the top
        if ~isempty(strfind(st.structure_id_path{ann}, isoCtxId))
            % if the track didn't get to cortex yet, keep looking...
            gotToCtx = true;
        end
    end
    ProbePoints{end+1}=[start_point(1), start_point(3), start_point(2)]; % this is actually the point where we exit the brain
    % plot brain grid
    if ~any(ishandle(brainfig))
        fwireframe = plotBrainGrid([], [], fwireframe, black_brain); hold on;
        fwireframe.InvertHardcopy = 'off';
        brainfig = gcf;
    end
    % plot probe points
    figure(brainfig)
    %hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);
    
    % plot brain entry point
    %plot3(m(1), m(3), m(2), 'k*','linewidth',1)
    
    % use the deepest clicked point as the tip of the probe, if no scaling provided (scaling_factor = false)
    if use_tip_to_get_reference_probe_length
        % find length of probe in reference atlas space
        [depth, tip_index] = max(curr_probePoints(:,2));
        reference_probe_length_tip = sqrt(sum((curr_probePoints(tip_index,:) - start_point).^2));
        
        % and the corresponding scaling factor
        shrinkage_factor = (reference_probe_length_tip / 100) / probe_length;
        
        % display the scaling
        disp(['probe length of ' num2str(reference_probe_length_tip/100) ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
        disp(['probe scaling of ' num2str(shrinkage_factor)]); disp(' ');
        
        % plot line the length of the probe in reference space
        probe_length = round(reference_probe_length_tip);
        
        % if scaling_factor is user-defined as some numer, use it to plot the length of the probe
    else
        probe_length = round(reference_probe_length * 100);
    end
    
    % find the percent of the probe occupied by electrodes
    percent_of_tract_with_active_sites = min([active_probe_length / probe_length, 1.0]);
    active_site_start = probe_length*(1-percent_of_tract_with_active_sites);
    active_probe_position = round([active_site_start  probe_length]);
    %probe_length = probe_length*.24;
    % plot line the length of the active probe sites in reference space
    plot3(start_point(1)+p(1)*[active_probe_position(1) active_probe_position(2)], start_point(3)+p(3)*[active_probe_position(1) active_probe_position(2)], start_point(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 2);
    % plot line the length of the entire probe in reference space
    plot3(start_point(1)+p(1)*[1 probe_length], start_point(3)+p(3)*[1 probe_length], start_point(2)+p(2)*[1 probe_length], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 2, 'LineStyle',':');
    x=start_point(1)+p(1)*[active_probe_position(1) active_probe_position(2)];
    z=start_point(2)+p(2)*[active_probe_position(1) active_probe_position(2)];
    y=start_point(3)+p(3)*[active_probe_position(1) active_probe_position(2)];
    endPoints = [x(2),y(2),z(2)];
    ProbePoints{end}=cat(1,ProbePoints{end},endPoints);
end
end