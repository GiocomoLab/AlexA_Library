% filenames={'npF2_1015_contrasttrack_gainchanges_2.mat',...
%     'npF2_1016_contrasttrack_gainchanges_1.mat',...
%     'npF3_1018_contrasttrack_gainchanges_1.mat',...
%     'npF3_1019_contrasttrack_gainchanges_contrast_1.mat',...
%     'npF4_1023_gaincontrast_1.mat',...
%     'npF4_1025_gaincontrast_2.mat '};
%restoredefaultpath
if ispc()
addpath(genpath('C:\code\AlexA_Library'));
addpath(genpath('C:\code\boundedline'));
addpath(genpath('F:\code\cortexlab_spikes'));

% filenames = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
%     'G2/1211_mismatch_1/1211_mismatch_1.mat',...
%     'G2/1212_mismatch_1/1212_mismatch_1.mat',...
%     'G5/1207_mismatch_1/1207_mismatch_1.mat',...
%     'G5/1210_mismatch_1/1210_mismatch_1.mat'
%     };

filenames = dir('Z:\giocomo\attialex\NP_DATA\mismatch\*mismatch*.mat');
filenames = dir('Z:\giocomo\attialex\NP_DATA\*mismatch*.mat');

root_dir='F:\';
else
run('/home/users/attialex/AlexA_Library/default_paths.m')
filenames=dir(fullfile(OAK,'attialex','NP_DATA','*mismatch*.mat'));
end

beh_varlist={'AID_B','MMRun','RunOFF','MMAllRun'};
varlist={'AID','CGS','DEPTH','avgMM','avgRunOn','avgRunOff','CID','session_name','session_type','REGION','PARENT'};
%%
aggregateData=struct();
aggregateBeh=struct();
for ii =1:length(varlist)
    aggregateData.(varlist{ii}) = [];
    eval([varlist{ii} '= [];'])
end
for ii=1:length(beh_varlist)
    aggregateBeh.(beh_varlist{ii})={};
end
MM_snps={};
%%
session_table = readtable('Z:\giocomo\attialex\NP_DATA\data_summary_June2019.xlsx');
%session_table = readtable(fullfile(OAK,'attialex','NP_DATA','data_summary_June2019.xlsx'));
session_names = session_table.SessionName;
%idx = strcmp(filenames(1).name(1:end-4),session_names);
for iF=1:numel(filenames)
    try
    %clear all
    clear anatomy
    load(fullfile(filenames(iF).folder, filenames(iF).name));
    session_name{iF} = filenames(iF).name;
    idx = strcmp(filenames(1).name(1:end-4),session_names);
    if nnz(idx)==0
        session_type{iF}='';
    else
    session_type{iF} = session_table.SessionType{idx};
    end
    extractMismatchTimes;
    aggregateBeh.MMRun{iF}=(squeeze(adata));
    aggregateBeh.MMAllRun{iF}=squeeze(adata_all);
    aggregateBeh.RunOff{iF} =squeeze(adataROFF(1,:,:));
%     resp = squeeze(mean(mm_rate,2));
%     tmpRON=squeeze(mean(rON_rate,2));
%     tmpRON=tmpRON(:,1:20:end);
%     tmpROFF=squeeze(mean(rOFF_rate,2));
%     tmpROFF=tmpROFF(:,1:20:end);
 
    avgMM = cat(1,avgMM,resp(:,1:20:end));
    avgRunOn = cat(1,avgRunOn,tmpRON);
    avgRunOff = cat(1,avgRunOff,tmpROFF);
    AID = cat(1,AID,ones(length(sp.cgs),1)*iF);
    CGS = cat(1,CGS,sp.cgs');
    CID = cat(1,CID,sp.cids');
    [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
        templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
    depth=zeros(length(sp.cgs),1);
    for iC=1:length(sp.cgs)
        depth(iC)=mean(spikeDepths(sp.clu==sp.cids(iC)));
    end
    DEPTH = cat(1,DEPTH,depth);
    drawnow;
    %sprintf('Now working on: %s',filenames{iF})
    %corrMbyDepth
    %plot_pause_sequence
    if exist('anatomy','var')
        if isfield(anatomy,'region_shifted')
            tmp_region = anatomy.region_shifted;
            tmp_parent = anatomy.parent_shifted;
        else
            if isfield(anatomy,'cluster_region')
            tmp_region = anatomy.cluster_region;
            else % because for now we only have cluster parent for mec data
                tmp_region = anatomy.cluster_parent;
            end
            
            tmp_parent = anatomy.cluster_parent;
        end
          if numel(tmp_region) ~= numel(depth)
              pp=anatomy.cluster_parent(sp.cids);
              tmp_parent = pp;
              tmp_
              error('anatomy and real clusters do not match')
          end
    else
        tmp_region = cell(1,numel(depth));
        tmp_parent = cell(1,numel(depth));
    end
    if ~isrow(tmp_parent)
        tmp_parent = tmp_parent';
    end
    if ~isrow(tmp_region)
        tmp_region = tmp_region';
    end
    PARENT=cat(2,PARENT,tmp_parent);
    REGION = cat(2,REGION,tmp_region);

catch ME
    warning(sprintf('error with File %s',filenames(iF).name));
    display(ME.message)
    end
    end
for ii =1:length(varlist)
    aggregateData.(varlist{ii}) = eval(varlist{ii});
end
aggregateData.MM_snps=MM_snps;
% aggregateData.avgMM=avgMM;
% aggregateData.avgRun=
% aggregateData.AID=AID;
% aggregateData.CGS=CGS;
% aggregateData.DEPTH = DEPTH;
%clearvars -except agg* avgMM filenames root_dir AID
save(fullfile(OAK,'attialex',strcat('MM_aggregate',date,'.mat')),'aggregateData','-v7.3')
aggregateData.MM_snps=[];
save(fullfile(OAK,'attialex',strcat('MM_aggregateSmall',date,'.mat')),'aggregateData','-v7.3')

%%
params=struct();
params.winIDX=-200:200;
params.masterTime=params.winIDX/50;
params.xLim=[-1 3];
figure

plotAVGSEM(aggregateData.avgMM(aggregateData.CGS>1,:)',gca,'parameters',params,'ms',true,'baseline',165:190)
xlabel('time [s]')
ylabel('Firing rate change [Hz]')
grid on

[a,b]=sort(aggregateData.DEPTH);
figure
plotAVGSEM(aggregateData.avgMM(b(1:200),:)',gca,'parameters',params,'ms',true,'baseline',165:190)
plotAVGSEM(aggregateData.avgMM(b(end-200:end),:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])
xlabel('time [s]')
ylabel('Firing rate change [Hz]')
grid on
legend({'ventral','Dorsal'})


%%
figure
plotAVGSEM(aggregateData.avgRunOff(aggregateData.CGS==2,:)',gca,'parameters',params,'ms',true,'baseline',165:190)
hold on
plotAVGSEM(aggregateData.avgMM(aggregateData.CGS==2,:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])
legend({'Run Off','MM'})
xlabel('time [s]')
ylabel('Firing rate change [Hz]')
grid on
%figure


%%
% [a,~]=max(aggregateData.avgMM(aggregateData.CGS==2,150:350),[],2);
% ff_plot=bsxfun(@rdivide,aggregateData.avgMM(aggregateData.CGS==2,150:350),a);
% figure
% imagesc(ff_plot)
% set(gca,'XTick',[0:50:200])
% set(gca,'XTickLabel',[-1:1:4])
% colormap summer
% colorbar
%set(gca,'
%%
[a,b]=sort(aggregateData.DEPTH(aggregateData.CGS==2));
%ff_plot=bsxfun(@rdivide,aggregateData.avgMM(aggregateData.CGS==2,150:350),a);
ff_plot=aggregateData.avgMM-mean(aggregateData.avgMM(:,170:195),2);
ff_plot=ff_plot(aggregateData.CGS==2,:);
figure
imagesc(ff_plot(:,150:350),[-10 10])
set(gca,'XTick',[0:50:200])
set(gca,'XTickLabel',[-1:1:4])
colormap jet
colorbar
xlabel('Time, [s]')
ylabel('Unit #')
%%
[a,b]=sort(aggregateData.DEPTH(aggregateData.CGS==2));
%ff_plot=bsxfun(@rdivide,aggregateData.avgMM(aggregateData.CGS==2,150:350),a);
ff_plot=aggregateData.avgMM-mean(aggregateData.avgMM(:,170:195),2);
ff_plot=ff_plot(aggregateData.CGS==2,:);
figure
imagesc(ff_plot(b,150:350),[-10 10])
set(gca,'XTick',[0:50:200])
set(gca,'XTickLabel',[-1:1:4])
set(gca,'YTick',[0:100:1000],'YTickLabel',[3000:-300:0])
colormap jet
colorbar
xlabel('Time, [s]')
ylabel('Distance from Tip')
%ylabel('a')
%%
CID=aggregateData.CGS==2 & ~isnan(aggregateData.Stability);
[a,b]=sort(aggregateData.Stability(CID),'desc');

tmpMM=aggregateData.avgMM(CID,:);

figure

plotAVGSEM(tmpMM(b(1:100),:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])
plotAVGSEM(tmpMM(b(end-100:end),:)',gca,'parameters',params,'ms',true,'baseline',165:190)

legend({'Low Stabilitly','High Stability'})
xlabel('time [s]')
ylabel('Firing rate change [Hz]')
grid on
%%
figure
legs={};
run_mat = [];
mm_mat = [];
for ii=1:5
    %size(aggregateBeh.RunOff{ii})
    plot([-4:0.02:4],mean(aggregateBeh.MMRun{ii}))
    run_mat=cat(1,run_mat,mean(aggregateBeh.MMRun{ii}));
    IDX=aggregateData.AID==ii & aggregateData.CGS==2;
    mm_mat=cat(1,mm_mat,mean(aggregateData.avgMM(IDX,:)));
    hold on
    tmp = strsplit(filenames{ii},'/');
    legs{ii}=tmp{1};
end
legend(legs)
grid on
ylabel('run speed')
xlabel('time')
figure
subplot(1,2,1)
imagesc(run_mat)
set(gca,'XTick',[0:50:400])
set(gca,'XTickLabel',[-4:1:4])
title('MM Behavior')
colorbar
xlabel('time')
subplot(1,2,2)
imagesc(mm_mat)
set(gca,'XTick',[0:50:400])
set(gca,'XTickLabel',[-4:1:4])
title('MM Response')
colorbar
xlabel('time')

%%
figure
legs={};
run_mat = [];
mm_mat = [];
for ii=1:5
    %size(aggregateBeh.RunOff{ii})
    plot(mean(aggregateBeh.MMRun{ii}))
    run_mat=cat(1,run_mat,mean(aggregateBeh.MMRun{ii}));
    IDX=aggregateData.AID==ii & aggregateData.CGS==2;
    mm_mat=cat(1,mm_mat,mean(aggregateData.avgMM(IDX,:)));
    hold on
    tmp = strsplit(filenames{ii},'/');
    legs{ii}=tmp{1};
end
legend(legs)
figure
subplot(1,2,1)
imagesc(run_mat)
subplot(1,2,2)
imagesc(mm_mat)
%%
figure
legs={};
run_mat = [];
mm_mat = [];
colors = winter(5);
for ii=1:5
    %figure
    IDX=aggregateData.AID==ii & aggregateData.CGS>=1;
    nnz(IDX)
    plotAVGSEM(aggregateData.avgMM(IDX,:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',colors(ii,:))
end
ylabel('Change in Firing Rate [Hz]')
xlabel('Time [s]')
