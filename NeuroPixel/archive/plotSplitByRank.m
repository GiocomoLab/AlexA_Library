function plotSplitByRank(spikes,runTrace,groupLabel,nSplit)

%LLI_variables;

ranks=[];
[a]=unique(groupLabel);
runAll=[];
actAll=[];
a=a';

% for id=a
%     idx=id==groupLabel;
%     runs=-1*cat(2,auxCell{idx});
%     act=cat(2,actCell{idx});
%     runAll=cat(2,runAll,runs);
%     actAll=cat(2,actAll,act);
%     
%     runS=mean(runs(85:115,:));
%     [~,~,cOff]=unique(runS);
%     rankSite=cOff/max(cOff)';
%     ranks=cat(1,ranks,rankSite);
% end



frac=1/nSplit;

all=figure;

col=jet(nSplit);
for ii=1:nSplit
    lims=frac*[ii-1 ii]
    tmpidx=ranks>lims(1) & ranks<=lims(2);
    
    [m1,sem1]=getMSEM(runAll(:,tmpidx));
    [m2,sem2]=getMSEM(actAll(:,tmpidx));
    

        m2=m2-mean(m2(lli.baseline));
 
    
    %     figure(single)
    %     subplot(nsplits,2,(ii-1)*2+1)
    %     boundedline(t1,m1,sem1)
    %     subplot(nsplits,2,ii*2)
    %     boundedline(t2,m2,sem2);
    
    figure(all)
    subplot(121)
    hold on
    boundedline(lli.masterTime,m1,sem1,'cmap',col(ii,:),'transparency',.2);
    %plot(t1,m1,'Color',col(ii,:))
    subplot(122)
    hold on
    %plot(t2,m2,'Color',col(ii,:));
    boundedline(lli.masterTime,m2,sem2,'cmap',col(ii,:),'transparency',.2);
    
end
end

    function [m,sem]=getMSEM(trace)
        m=nanmean(trace,2);
        sem=nanstd(trace,[],2)/sqrt(size(trace,2));
    end


