data_dir = fullfile('/oak/stanford/groups/giocomo/attialex','NP_DATA');
session_name = dir(fullfile(data_dir,'AA*gain*.mat'));
nS=numel(session_name);
AVG_SPEED = zeros(nS,401);
%
for iS=1:numel(session_name)
%load session
load(fullfile(session_name(iS).folder,session_name(iS).name))
%calc speed
params=readtable('UniversalParams.xlsx');
speed = calcSpeed(posx,params);

%% trial_idx (baseline only)
bl_trial = find(trial_gain==1 & trial_contrast==100);
bl_trial(bl_trial>max(trial))=[];
valid_idx = ismember(trial,bl_trial);
posx(posx<0)=0;
posx(posx>400)=400;
posx_discrete = discretize(posx,0:401);
speed_mat = zeros(max(trial),401);
for iT=1:max(trial)
    for iPos = 1:401
        idx = trial == iT & posx_discrete == iPos;
        tmp = mean(speed(idx));
        speed_mat(iT,iPos)=tmp;
    end
end
avg_speed = nanmean(speed_mat(bl_trial,:));
% for iP = 1:401
%     idx = posx_discrete == iP & valid_idx;
%     avg_speed(iP)=mean(speed(idx));
% end
AVG_SPEED(iS,:)=avg_speed;
end
%%
params=struct();
params.winIDX=1:401;
params.masterTime=0:400;
params.xLim=[0 400];
set(0,'DefaultFigureRenderer','painters')
figure('Position',[1640         764         560         214]);
plotAVGSEM(AVG_SPEED',gca,'parameters',params,'ms',false)
saveas(gcf,'/oak/stanford/groups/giocomo/attialex/behav.pdf')