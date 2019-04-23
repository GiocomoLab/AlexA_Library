

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
        depth=zeros(length(sp.cgs),1);
        for iC=1:length(sp.cgs)
            depth(iC)=median(spikeDepths(sp.clu==sp.cids(iC)));
        end
        %%
        [a,sid]=sort(depth);
        a=a-a(1);
        superficial=sp.cgs;
        idx=ismember(sp.clu,[1 3]);
        cluID=zeros(size(sp.clu));
        nsplits=12;
        stepsize=a(end)/nsplits;
        
        for isplit=1:nsplits
            start = (isplit-1)*stepsize;
            stop = isplit*stepsize;
            current = sp.cids(sid(a>=start & a<=stop));
            tmpidx=ismember(sp.clu,current);
            cluID(tmpidx)=isplit;
            
        end
            
        %%
        [CGR,b]=CCG(sp.st,cluID,'binSize',[0.001],'duration',[0.2]);
        
        %%
        mean_depth = zeros(1,nsplits);
        for ii=1:nsplits
            mean_depth(ii)=mean(spikeDepths(cluID==ii));
        end

        %%
        figure
        hold on
        leg={};
        delays = zeros(1,nsplits);
        for ii=2:nsplits
            tmp=squeeze(CGR(:,1,ii));
            tmp=tmp/max(tmp);
            [~,maxi]=max(tmp);
            delays(ii)=b(maxi);
            plot(b,tmp)
            leg{end+1}=sprintf('%.2f',mean_depth(ii));
        end
        legend(leg)
        
        %%
       
        