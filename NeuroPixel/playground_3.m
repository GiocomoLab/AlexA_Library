
data_dir = 'C:\Users\giocomolab\Dropbox\Work\neuropixels\data\';
session_name = {'npF2_1015_contrasttrack_gainchanges_2', ...
    'npF2_1016_contrasttrack_gainchanges_1', ...
    'npF3_1018_contrasttrack_gainchanges_1',...
    'npF3_1019_contrasttrack_gainchanges_contrast_1',...
    'npF4_1023_gaincontrast_1', ...
    'npF4_1025_gaincontrast_2',...
    'npG1_1218_gaincontrast_2',...
    'npG1_1219_gain_1',...
    'npG2_1211_gain_1',...
    'npG2_1212_gaincontrast_1',...,
    'npG2_1213_gain_1',...
    'npG2_1214_gaincontrast_1',...
    'npG4_1203_gain_2',...
    'npG4_1204_gaincontrast_2',...
    'npG5_1206_gaincontrast_1',...
    'npG5_1207_gain_1',...
    'npG5_1210_gaincontrast_2'};

session_num=7;
load(fullfile(data_dir,strcat(session_name{session_num},'.mat')));

trials=[1:max(trial)];
spC=[];
dwell_time=[];
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    edges=[0:2:402];
    edges(1)=-15;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges,'max_clust',max(sp.cids));
    spC=cat(3,spC,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
spC=spC(sp.cids+1,:,:);
dt=dwell_time';

dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spC,1)
    spC(ii,:,:)=spC(ii,:,:)./dt;
end
spC=spC/0.02;

%%
% sp.posx=posx;
% sp.trial=trial;
% sp.post=post;
% psthViewer2_Alex(sp, eventTimes, window, trialGroups,spC);


%%
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(fullfile('F:',session_name{session_num}(3:4),session_name{session_num}(6:end)));
%%
mkdir(fullfile('F:','tmp_images',session_name{session_num},'raster_amplitude'))

x_width=12.25 ;y_width=8.125;
ff=figure('PaperUnits', 'inches','PaperPosition', [0 0 x_width y_width]);    
plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
xlim([min(sp.st), max(sp.st)])
saveas(ff,fullfile('F:','tmp_images',session_name{session_num},['drift_map_' session_name{session_num} '.png']),'png');

%%

x_width=12.25 ;y_width=8.125;
ff=figure('PaperUnits', 'inches','PaperPosition', [0 0 x_width y_width]);    

for ii =1:length(sp.cids)
    if sp.cgs(ii)==2
idx = sp.cids(ii);
h=subplot(2,1,1);
plot(sp.st(sp.clu==idx),sp.tempScalingAmps(sp.clu==idx),'.','MarkerSize',8)
title([session_name{session_num} ' cluster ' num2str(sp.cids(ii))],'Interpreter', 'none');
xlim([min(sp.st), max(sp.st)])
hold off;
xlabel('time [s]')
ylabel('Ampl [au]')

box off;

h=subplot(2,1,2);
hold off;
fr=squeeze(spC(ii,:,:));
imagesc(fr)

set(h,'Clipping','Off')
axes(h)
for iT=20:20:200
    tt=min(post(trial==iT));
    frac = tt/max(post);
    trialX=round(frac*max(trial));
    line([iT trialX],[0 -80],'Color','red','LineStyle','--');
end
xlabel('trial #')
ylabel('position')

saveas(ff,fullfile('F:','tmp_images',session_name{session_num},'raster_amplitude',...
            [session_name{session_num} sprintf('_cell_%d.png',sp.cids(ii))]),'png');
        %clf(ff);
    end
end