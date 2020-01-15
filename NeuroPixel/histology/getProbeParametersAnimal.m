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
elseif strcmp(animalName,'AA_190830_044')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190830_044\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.8 3.4 3.65 3.7]; 
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190830_044\combined\processed';
elseif strcmp(animalName,'AA_190830_045')
    image_save_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190830_045\combined';
    probe_save_name_suffix = 'combined';
    probe_lengths = [2 3.2 3.0 3.7 3.8]; 
    processed_images_folder = 'Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190830_045\combined\processed';
elseif strcmp(animalName,'AA_190830_047')
    image_save_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined');
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.75 3.9 3.8 3.7]; 
    processed_images_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined\processed');
elseif strcmp(animalName,'AA_190830_046')
    image_save_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined');
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.7 3.9 3.7 3.9 3.8]; 
    processed_images_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined\processed');

elseif strcmp(animalName,'AA_190906_049')
    image_save_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined');
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.7]; 
    processed_images_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined\processed');

elseif strcmp(animalName,'AA_190906_050')
    image_save_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined');
    probe_save_name_suffix = 'combined';
    probe_lengths = [3.7 3.7 3.7 3.9 3.8]; 
    processed_images_folder = strcat('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\',animalName,'\combined\processed');
end

end