function [filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,data_dir)

files = dir(fullfile(data_dir,'*.mat'));
filenames={};
triggers = {};
for iF = 1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name),'anatomy','trial_gain','trial_contrast');
    gaincontrastcombo = false;
    if contains(files(iF).name,'mismatch') || contains(files(iF).name,'playback') || contains(files(iF).name,'dark')
        continue
    end
    if isfield(data,'trial_gain') && isfield(data,'trial_contrast')
        if strcmp(gain_to_look_at,'baseline') || any(gain_to_look_at == 1)
            %find condition with given contrast and no gain change
            bin = data.trial_gain == gain_to_look_at & data.trial_contrast == contrast;
            trigger_tmp = strfind(bin',ones(1,16)) +8;
            trigger_tmp(trigger_tmp <10) = [];
            gaincontrastcombo = nnz(trigger_tmp)>0;
        elseif any(gain_to_look_at == 0)
            %special case: find where gain ==1 and for a given contrast
            %change
            gaincontrastcombo = any(data.trial_gain == 1 & data.trial_contrast == contrast);
            trigger_tmp = strfind((data.trial_gain == 1 & data.trial_contrast == contrast)',[0 1])+1;
        elseif isempty(gain_to_look_at)
            %find contrast steps (100 to contrast)
            pot = strfind(data.trial_contrast'  == contrast,[0 1])+1;
            pot = pot(data.trial_contrast(pot-1)==100);
            gaincontrastcombo = any(pot);
            trigger_tmp = pot;
         
        else
        gaincontrastcombo = any(data.trial_gain == gain_to_look_at & data.trial_contrast == contrast);
        trigger_tmp = strfind((data.trial_gain == gain_to_look_at & data.trial_contrast == contrast)',[0 1])+1;
        end
        end
    containsregion = false;
    if isfield(data,'anatomy')
        containsregion = any(startsWith(data.anatomy.cluster_parent,region));
    end
    
    if gaincontrastcombo && containsregion

        filenames{end+1}=fullfile(files(iF).folder,files(iF).name);
        triggers{end+1}=trigger_tmp;
    end
end
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

end

