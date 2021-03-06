root = 'Z:\giocomo\attialex\NP_DATA\';
matfiles = dir(fullfile(root,'npI1_0417_dark_1.mat'));

pairs = [];
AID = [];
depth_pre=[];
depth_post=[];
CGRs=[];
for iF=1:length(matfiles)
    clear connected
    load(fullfile(root,matfiles(iF).name),'connected')
    
    if exist('connected','var')
        vars = {'sp','PCausal'};
        load(fullfile(root,matfiles(iF).name),vars{:})
        good_cells = sp.cids(sp.cgs==2);
        
        [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
        depth=zeros(length(sp.cgs),1);
        for iC=1:length(sp.cgs)
            depth(iC)=median(spikeDepths(sp.clu==sp.cids(iC)));
        end
        
        idx=ismember(sp.clu,connected(:));
        [CGR,b]=CCG(sp.st(idx),double(sp.clu(idx))+1,'binSize',[0.0004],'duration',[0.2]);
        pairs = cat(1,pairs,connected);
        for ii=1:length(connected)
            if PCausal(connected(ii,1)+1,connected(ii,2)+1)<0.05
            depth_pre(end+1)=depth(good_cells==connected(ii,1));
            depth_post(end+1) = depth(good_cells==connected(ii,2));
            AID(end+1)=iF;
            CGRs=cat(2,CGRs,squeeze(CGR(:,connected(ii,2)+1,connected(ii,1)+1)));
            end
        end
    else
        sprintf('%s does not contain connected',matfiles(iF).name)
    end
    
    
end
%%
figure
plot([depth_pre',depth_post']')
set(gca,'XTick',[1 2],'XTickLabel',{'Pre', 'Post'})
xlim([.9 2.1])
ylabel('distance from tip')

%%
distance = abs(depth_pre-depth_post);

slices =10;
[a,id]=sort(distance);
stepsize=floor(length(id)/slices);
figure
dd={};
for iS=1:slices
    idx=id((iS-1)*stepsize+1:iS*stepsize);
    slice=CGRs(:,idx);
    mm=max(slice,[],1);
    tmp=slice./mm;
    avg=mean(tmp,2);
    hold on
    plot([-.1:0.0004:.1],avg)
    dd{iS}=sprintf('%d',mean(distance(idx)));
end
legend(dd)
%%
