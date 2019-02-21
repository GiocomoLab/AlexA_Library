%% normalize by baseline
session_name = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
    'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    'G5/1210_mismatch_1/1210_mismatch_1.mat'
    };


root_dir='F:\';

origin = nan(numel(session_name),3); % for parameterizing probe track
unit_vector = nan(numel(session_name),3); % for parameterizing probe track
figure
hold on
for session_num = 1:numel(session_name)
    clear histology
    load(fullfile(root_dir,session_name{session_num}));
    if exist('histology','var')
    vec1 = histology.MEC_entry-histology.probe_term;
    vec1(3) = -vec1(3); % flip Z coord
    vec1 = vec1./repmat(sqrt(sum(vec1.^2,2)),1,3); % normalize
    vec2 = histology.probe_term; 
    vec2(3) = -vec2(3); % flip Z coord
    unit_vector(session_num,:) = vec1;
    origin(session_num,:) = vec2;
    end

    %stability = stability_all{i};
    IDX=aggregateData.AID==session_num & aggregateData.CGS==2;
    spike_depth = aggregateData.DEPTH(IDX);
    MM_resp = mean(aggregateData.avgMM(IDX,205:255),2)./mean(aggregateData.avgMM(IDX,5:195),2);
    
    MM_resp=mat2gray(MM_resp)*10;
    
    % get 3D anatomical locations of units for this session
    pos3D = origin(session_num,:)+spike_depth.*unit_vector(session_num,:);
    
    % make 3D plot (maybe there's a better way than a for loop...)
    for j = 1:numel(spike_depth)
        if MM_resp(j)>0
        plot3(pos3D(j,1),pos3D(j,2),pos3D(j,3),'ko',...
            'MarkerSize',MM_resp(j));
        end
    end
end
%%
session_name = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
    'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    'G5/1210_mismatch_1/1210_mismatch_1.mat'
    };


root_dir='F:\';

origin = nan(numel(session_name),3); % for parameterizing probe track
unit_vector = nan(numel(session_name),3); % for parameterizing probe track
figure
hold on
for session_num = 1:numel(session_name)
    clear histology
    load(fullfile(root_dir,session_name{session_num}));
    if exist('histology','var')
    vec1 = histology.MEC_entry-histology.probe_term;
    vec1(3) = -vec1(3); % flip Z coord
    vec1 = vec1./repmat(sqrt(sum(vec1.^2,2)),1,3); % normalize
    vec2 = histology.probe_term; 
    vec2(3) = -vec2(3); % flip Z coord
    unit_vector(session_num,:) = vec1;
    origin(session_num,:) = vec2;
    end

    %stability = stability_all{i};
    IDX=aggregateData.AID==session_num & aggregateData.CGS==2;
    spike_depth = aggregateData.DEPTH(IDX);
    MM_resp = mean(aggregateData.avgMM(IDX,205:255),2)-mean(aggregateData.avgMM(IDX,145:195),2);
    MM_resp(MM_resp>5)=5;
    MM_resp(MM_resp<-5)=-5;
    MM_resp_size=mat2gray(abs(MM_resp))*20;
    
    % get 3D anatomical locations of units for this session
    pos3D = origin(session_num,:)+spike_depth.*unit_vector(session_num,:);
    
    % make 3D plot (maybe there's a better way than a for loop...)
%     for j = 1:numel(spike_depth)
%         if MM_resp(j)>0
%         plot3(pos3D(j,1),pos3D(j,2),pos3D(j,3),'o',...
%             'MarkerSize',MM_resp_size(j),'MarkerColor');
%         end
%     end
    iidx=MM_resp_size>0;

    scatter3(pos3D(iidx,1),pos3D(iidx,2),pos3D(iidx,3),MM_resp_size(iidx),MM_resp(iidx))
end
[x y] = meshgrid(-1000:50:1000);
z=zeros(size(x));
surf(x,y,z,ones(size(x,1),size(x,2),3)*.9,'FaceAlpha',0.5,'EdgeColor','none')
xlim([-500 500])
ylim([0 800])
grid on

%% extend aggregateData
session_name = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
    'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    'G5/1210_mismatch_1/1210_mismatch_1.mat'
    };
aggregateData.Hist_Depth = [];

root_dir='F:\';

origin = nan(numel(session_name),3); % for parameterizing probe track
unit_vector = nan(numel(session_name),3); % for parameterizing probe track
figure
hold on
for session_num = 1:numel(session_name)
    clear histology
    load(fullfile(root_dir,session_name{session_num}));
    if exist('histology','var')
    vec1 = histology.MEC_entry-histology.probe_term;
    vec1(3) = -vec1(3); % flip Z coord
    vec1 = vec1./repmat(sqrt(sum(vec1.^2,2)),1,3); % normalize
    vec2 = histology.probe_term; 
    vec2(3) = -vec2(3); % flip Z coord
    unit_vector(session_num,:) = vec1;
    origin(session_num,:) = vec2;
    end

    
    IDX=aggregateData.AID==session_num;
    spike_depth = aggregateData.DEPTH(IDX);
    
    
    
    pos3D = origin(session_num,:)+spike_depth.*unit_vector(session_num,:);
    aggregateData.Hist_Depth = cat(1,aggregateData.Hist_Depth,pos3D(:,3));
  
end
