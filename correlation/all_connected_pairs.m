root = 'Y:\giocomo\attialex\NP_DATA\';
matfiles = dir([root '*.mat']);

pairs = [];
AID = [];
depth_pre=[];
depth_post=[];
CGRs=[];
for iF=1:length(matfiles)
    clear connected
    load(fullfile(root,matfiles(iF).name),'connected')
    
    if exist('connected','var')
        load(fullfile(root,matfiles(iF).name),'sp')
        good_cells = sp.cids(sp.cgs==2);
        
        [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
        depth=zeros(length(sp.cgs),1);
        for iC=1:length(sp.cgs)
            depth(iC)=median(spikeDepths(sp.clu==sp.cids(iC)));
        end
        
        idx=ismember(sp.clu,connected(:));
        [CGR,b]=CCG(sp.st(idx),double(sp.clu(idx)),'binSize',[0.0004],'duration',[0.2]);
        pairs = cat(1,pairs,connected);
        for ii=1:length(connected)
            depth_pre(end+1)=depth(good_cells==connected(ii,1));
            depth_post(end+1) = depth(good_cells==connected(ii,2));
            AID(end+1)=iF;
            CGRs=cat(2,CGRs,squeeze(CGR(:,connected(ii,1),connected(ii,2))));
            
        end
    end
    
end
%%
figure
plot([depth_pre',depth_post']')

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