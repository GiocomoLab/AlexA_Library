matfiles = dir('Z:\giocomo\export\data\Projects\AlexA_NP\AA_2009*\ks_data\*\AA*.mat');

%%
for iF=1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    [tmp]=fileparts(matfiles(iF).folder);
    vid_folder = fileparts(tmp);
    vid_file = fullfile(vid_folder,[sn '_proc.mat']);
    if ~isfile(vid_file)
        fprintf('no video file for %s \n',sn)
        continue
    end
    ft_file = fullfile(vid_folder,[sn '.ft']);
    ft_data = importdata(ft_file);
    vid_data = load(vid_file);
    tmp_file = strrep(vid_file,'_proc','_proc_area');
    if isfile(tmp_file)
        p=load(tmp_file);
        area_raw = p.area;
    else
        fprintf('no raw area for %s \n',sn)
        area_raw = nan(size(vid_data.pupil{1}.area));
    end
    pupil_area = vid_data.pupil{1}.area;
    pupil_com = vid_data.pupil{1}.com;
    pupil_area_raw = area_raw;
    vid_length = ft_data(end,2)-ft_data(1,2);
    
    %data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    %data = matfile(fullfile(matfiles(iF).folder,matfiles(iF).name));
    data = load(fullfile(matfiles(iF).folder,matfiles(iF).name),'post');
    rec_length = data.post(end)-data.post(1);
    time_diff = vid_length-rec_length;
    if abs(time_diff)>0.2
        fprintf('large diff (%.2f) with %s \n',time_diff,sn)
        continue
    end
    
    vid_frame_times =  ft_data(:,2)-ft_data(1,2)+(-1*time_diff*0.5); %transforming vid_frame times by adding half of the difference between vr and video 
    pupil_area = interp1(vid_frame_times,pupil_area,data.post);
    pupil_com = interp1(vid_frame_times,pupil_com,data.post);
    pupil_area_raw = interp1(vid_frame_times,pupil_area_raw,data.post);
    save(fullfile(matfiles(iF).folder,matfiles(iF).name),'pupil_area','pupil_com','pupil_area_raw','-append')
    copyfile(fullfile(matfiles(iF).folder,matfiles(iF).name),'F:\Alex\matfiles_processed')
    %save(fullfile('F:\Alex\matfiles_new',matfiles(iF).name),'pupil_area','pupil_com','-append')

end
%%

