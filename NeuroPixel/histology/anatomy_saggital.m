%%
addpath(genpath('/home/users/attialex/AlexA_Library'));
addpath(genpath('/home/users/attialex/spikes/'));
%%
histology = process_histology();
session_summary = readtable('/oak/stanford/groups/giocomo/attialex/NP_DATA/session_summary.xlsx');
nSessions = size(session_summary,1);
data_path = '/oak/stanford/groups/giocomo/attialex/NP_DATA/mismatch';
tmp = datestr(histology.Date,'mmdd');
datestrings = {};
for ii=1:size(tmp,1)
    datestrings{ii}=tmp(ii,:);
end
%for each session
for ii=1:168
    try
        clear anatomy;
        clear dataset
        clear cluster_parent
        dataset = load(fullfile(data_path,session_summary.SessionName{ii}));
        
        [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(dataset.sp.temps, dataset.sp.winv, dataset.sp.ycoords, dataset.sp.spikeTemplates, dataset.sp.tempScalingAmps);
        depth=zeros(length(dataset.sp.cgs),1);
        for iC=1:length(dataset.sp.cgs)
            depth(iC)=mean(spikeDepths(dataset.sp.clu==dataset.sp.cids(iC)));
        end
        
        animal = strsplit(session_summary.SessionName{ii},'_');
        edate = animal{2};
        animal = animal{1};
        aidx = strcmp(histology.animal,animal);
        didx = strcmp(datestrings,edate);
        pp=find(aidx & didx');
        %display(sprintf('Name: %s, main: %s',session_summary.SessionName{ii},[animal, '_' ,histology.Session{aidx & didx'}]))
        
        anatomy.tip_distance = depth;
        cluster_parent = cell(numel(depth),1);
        if ~isempty(pp)
            anatomy.z2=histology.Z2(pp);
            anatomy.FinalDepth = histology.FinalDepth(pp);
            if ~isnan(anatomy.z2)
                for ic=1:length(depth)
                    if depth(ic)<anatomy.z2
                        cluster_parent{ic} = 'MEC';
                    else
                        cluster_parent{ic} = 'ECT';
                    end
                end
            else
                disp('using surface')
                %keyboard
                for ic=1:length(depth)
                    if depth(ic)<anatomy.FinalDepth - 1000
                        cluster_parent{ic} = 'MEC';
                    else
                        cluster_parent{ic} = 'ECT';
                    end
                end
            end
            
            anatomy.cluster_parent = cluster_parent;

        end
        save(fullfile(data_path,session_summary.SessionName{ii}),'anatomy','-append');
    catch ME
        disp(ME.message)
    end
           if exist('cluster_parent','var') & numel(unique(cluster_parent))==1
                disp(['only 1 region for ' session_summary.SessionName(ii)])
                error('ieie')
            end
end

