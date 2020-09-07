    params=struct();
    params.winIDX=-4000:4000;
    params.masterTime=params.winIDX/1000;
    params.xLim=[-1 3];
    paramsAdata=struct();
    paramsAdata.winIDX=-200:200;
    paramsAdata.masterTime=paramsAdata.winIDX/50;
    paramsAdata.xLim=[-1 3];
    allMM={};
        nSplit=4;

    for iS=nSplit
        allMM{iS}=[];
        allBeh{iS}=[];
    end
for iA=1:5 
    runS=mean(aggregateBeh.MMAllRun{iA}(:,175:225),2);
    [~,~,cOff]=unique(runS);
    rank=cOff/max(cOff)';
    frac=1/nSplit;
    kk=reshape(gausswin(201),1,[]);
    figure('Name',filenames{iF});
    col=jet(nSplit);
    for ii=1:nSplit
        lims=frac*[ii-1 ii];
        tmpidx=rank>lims(1) & rank<=lims(2);
        IDX=aggregateData.CGS(aggregateData.AID==iA)==2;
        
        resp = squeeze(mean(aggregateData.MM_snps{iA}(IDX,tmpidx,:),2));
        resp = convn(resp,kk,'same');
        subplot(2,1,1)
        plotAVGSEM(resp',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',col(ii,:))
        subplot(2,1,2)
        hold on
        plot(-4:0.02:4,squeeze(mean(aggregateBeh.MMAllRun{iA}(tmpidx,:))),'Color',col(ii,:))
        allMM{ii}=cat(1,allMM{ii},resp);
        allBeh{ii}=cat(1,allBeh{ii},mean(aggregateBeh.MMAllRun{iA}(tmpidx,:)));
    end
end

%%
figure
hold on
for iS=1:nSplit
    subplot(2,1,1)
    hold on
    plotAVGSEM(allMM{iS}',gca,'parameters',params,'ms',true,'baseline',3500:3999,'col',col(iS,:))
    subplot(2,1,2)
    hold on
    plot([-4:0.02:4],mean(allBeh{iS}),'Color',col(iS,:))
    xlim(params.xLim)
    
end
    
    