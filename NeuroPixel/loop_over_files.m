% filenames={'npF2_1015_contrasttrack_gainchanges_2.mat',...
%     'npF2_1016_contrasttrack_gainchanges_1.mat',...
%     'npF3_1018_contrasttrack_gainchanges_1.mat',...
%     'npF3_1019_contrasttrack_gainchanges_contrast_1.mat',...
%     'npF4_1023_gaincontrast_1.mat',...
%     'npF4_1025_gaincontrast_2.mat '};

filenames = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
    'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    'G5/1210_mismatch_1/1210_mismatch_1.mat'
    };


root_dir='F:\';
aggregateData=struct();
avgMM=[];
avgRunOn=[];
avgRunOff=[];

AID=[];
CGS=[];
DEPTH=[];
varlist={'AID','CGS','DEPTH','avgMM','avgRunOn','avgRunOff'};
for iF=1:5
    %clear all
    load([root_dir filenames{iF}]);
    tic;plot_mismatch_sequence;toc;
    resp = squeeze(mean(mm_rate,2));
    tmpRON=squeeze(mean(rON_rate,2));
    tmpRON=tmpRON(:,1:20:end);
    tmpROFF=squeeze(mean(rOFF_rate,2));
    tmpROFF=tmpROFF(:,1:20:end);
    
    avgMM = cat(1,avgMM,resp(:,1:20:end));
    avgRunOn = cat(1,avgRunOn,tmpRON);
    avgRunOff = cat(1,avgRunOff,tmpROFF);
    AID = cat(1,AID,ones(length(sp.cgs),1)*iF);
    CGS = cat(1,CGS,sp.cgs');
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
end
for ii =1:length(varlist)
    aggregateData.(varlist{ii}) = eval(varlist{ii});
end
% aggregateData.avgMM=avgMM;
% aggregateData.avgRun=
% aggregateData.AID=AID;
% aggregateData.CGS=CGS;
% aggregateData.DEPTH = DEPTH;
%clearvars -except agg* avgMM filenames root_dir AID

%%
params=struct();
params.winIDX=-200:200;
params.masterTime=params.winIDX/50;
params.xLim=[-1 3];
figure

plotAVGSEM(aggregateData.avgMM(aggregateData.CGS==2,:)',gca,'parameters',params,'ms',true,'baseline',165:190)

[a,b]=sort(aggregateData.DEPTH);
figure
plotAVGSEM(aggregateData.avgMM(b(1:round(length(b)/2)),:)',gca,'parameters',params,'ms',true,'baseline',165:190)
plotAVGSEM(aggregateData.avgMM(b(round(length(b)/2)):end,:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])


%%
figure
plotAVGSEM(aggregateData.avgRunOff(aggregateData.CGS==2,:)',gca,'parameters',params,'ms',true,'baseline',165:190)
hold on
plotAVGSEM(aggregateData.avgMM(aggregateData.CGS==2,:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])
legend({'Run Off','MM'})

%%
[a,~]=max(aggregateData.avgMM,[],2);
ff_plot=bsxfun(@rdivide,aggregateData.avgMM,a);
figure
imagesc(ff_plot)
%%
CID=aggregateData.CGS==2 & ~isnan(aggregateData.Stability);
[a,b]=sort(aggregateData.Stability(CID),'desc');

tmpMM=aggregateData.avgMM(CID,:);

figure

plotAVGSEM(tmpMM(b(1:300),:)',gca,'parameters',params,'ms',true,'baseline',165:190,'col',[1 0 0])
plotAVGSEM(tmpMM(b(301:end),:)',gca,'parameters',params,'ms',true,'baseline',165:190)



