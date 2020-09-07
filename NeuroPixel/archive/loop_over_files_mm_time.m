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
varlist={'AID','CGS','avgMM','SIG_VAL','CID','session_name','REGION','PARENT'};
%%
aggregateData=struct();
for ii =1:length(varlist)
    aggregateData.(varlist{ii}) = [];
    eval([varlist{ii} '= [];'])
end
for ii=1:length(beh_varlist)
    aggregateBeh.(beh_varlist{ii})={};
end
MM_snps={};
%%
%session_table = readtable('Z:\giocomo\attialex\NP_DATA\data_summary_June2019.xlsx');
%session_names = session_table.SessionName;
%idx = strcmp(filenames(1).name(1:end-4),session_names);
for iF=12:numel(filenames)
    %clear all
    clear anatomy
    load(fullfile(filenames(iF).folder, filenames(iF).name));
    session_name{iF} = filenames(iF).name;
    %idx = strcmp(filenames(1).name(1:end-4),session_names);
    %session_type{iF} = session_table.SessionType{idx};
    mm_timing_perSession
    avgMM=cat(1,avgMM,squeeze(mean(rate_mat,2)));
    SIG_VAL = cat(1,SIG_VAL,sig_val);
    
    AID = cat(1,AID,ones(length(sp.cgs),1)*iF);
    CGS = cat(1,CGS,sp.cgs');
    CID = cat(1,CID,sp.cids');
    nCLU = nnz(sp.cgs>=1);
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
          if numel(tmp_region) ~= nCLU
              error('anatomy and real clusters do not match')
          end
    else
        tmp_region = cell(1,nCLU);
        tmp_parent = cell(1,nCLU);
    end
    if ~isrow(tmp_parent)
        tmp_parent = tmp_parent';
    end
    if ~isrow(tmp_region)
        tmp_region = tmp_region';
    end
    PARENT=cat(2,PARENT,tmp_parent);
    REGION = cat(2,REGION,tmp_region);
    
    
    drawnow;
    %sprintf('Now working on: %s',filenames{iF})
    %corrMbyDepth
    %plot_pause_sequence
end
%%
for ii =1:length(varlist)
    aggregateData.(varlist{ii}) = eval(varlist{ii});
end

save(fullfile(OAK,'attialex',strcat('MM_aggregate_timing',date,'.mat')),'aggregateData','-v7.3')

%%
RANKING = [];
figure
nA=unique(aggregateData.AID);
for iA=1:numel(nA)
IDX=aggregateData.AID==iA & strcmp(aggregateData.PARENT,'VISp')';
if nnz(IDX)==0
    continue
end
sig_val=aggregateData.SIG_VAL(IDX,:);
avgMM=aggregateData.avgMM(IDX,:);
mmresp = mean(avgMM(:,220:250),2)-mean(avgMM(:,155:185),2);
[~,sidx]=sort(mmresp,'descend');
r=1:length(sidx);
r(sidx)=r;
r=r/max(r);
RANKING=cat(1,RANKING,r');
crossings = nan(size(sig_val,1),1);
zero_bin = 201;
for iC= 1:size(avgMM,1)
    a=strfind(sig_val(iC,:)<0.05,[0 1 1])+1;
    a(a<zero_bin)=[];
    if ~isempty(a)
        crossings(iC)=a(1);
    end
end
frac = .2;
n=round(frac*size(avgMM,1));
subplot(4,4,iA)
tvec = linspace(-4,4,size(rate_mat,3));
%tmp = mean(squeeze(mean(rate_mat(sidx(1:n),:,:),2)));
tmp = mean(avgMM(sidx(1:n),:));
plot(tvec,tmp)
hold on
ff=round(nanmedian(crossings(sidx(1:n))));
%ff=round(ff)+win(1)-1;
plot(tvec(ff),tmp(ff),'r*')

tmp_c = crossings(sidx(1:n));
tmp_c(isnan(tmp_c))=[];
%tmp_c=tmp_c+win(1)-1;
plot(tvec(tmp_c),tmp(ff)*ones(1,numel(tmp_c)),'g.')
xlim([-.5 2])

end
%%
IDX=RANKING<0.2;

sig_val=aggregateData.SIG_VAL(IDX,:);
avgMM=aggregateData.avgMM(IDX,:);

crossings = nan(size(sig_val,1),1);
zero_bin = 201;
for iC= 1:size(avgMM,1)
    a=strfind(sig_val(iC,:)<0.05,[0 1 1])+1;
    a(a<zero_bin)=[];
    if ~isempty(a)
        crossings(iC)=a(1);
    end
end
n=round(frac*size(avgMM,1));
figure
tvec = linspace(-4,4,size(avgMM,2));
tmp = mean(avgMM,1);
plot(tvec,tmp)
hold on
ff=round(nanmedian(crossings));
%ff=round(ff)+win(1)-1;

tmp_c = crossings;
tmp_c(isnan(tmp_c))=[];
%tmp_c=tmp_c+win(1)-1;
plot(tvec(tmp_c),tmp(ff)*ones(1,numel(tmp_c)),'g.')
plot(tvec(ff),tmp(ff),'r*')
xlabel('Time from MM onset')
xlim([-.2 1.5])
grid on

figure
histogram(tvec(tmp_c),[0:0.02:2])
title('histogram of response times')
xlabel('time from MM onset')
grid on

%%
[N,edges] = histcounts(tvec(tmp_c),[0:0.02:2]);
centers = (edges(1:end-1)+edges(2:end))/2;
figure
plot(centers,N,'.')
[~,ii]=max(N);
centers(ii)
%%
x=zeros(1,length(tmp_c))';
figure
beeswarm(x,tvec(tmp_c)')
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

