function [image_save_folder,probe_save_name_suffix,probe_lengths,processed_images_folder]=getProbeParametersAnimal(animalName)

if strcmp(animalName,'AA_190709_1')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_1\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [3 3.5 3.5 3.5 3.0];
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_1\combined\processed';
elseif strcmp(animalName,'AA_190709_2')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_2\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.5 3.8 3.8 3.8 3.5 4.3];
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_2\combined\processed';
elseif strcmp(animalName,'AA_190709_3')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_3\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.0 3.2 3.1 3.8 3.7];
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_3\combined\processed';
elseif strcmp(animalName,'AA_190709_4')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_4\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.3 3.5 3.5 3.8 3.8 3.8 3.5];
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_4\combined\processed';
elseif strcmp(animalName,'AA_190709_5')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_5\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.3 3.5 3.5 3.8 3.8 3.8 3.5]; %only have real data for last two penetrations
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_5\combined\processed';
end
end