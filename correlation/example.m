good_cells = sp.cids(sp.cgs==2);
idx=ismember(sp.clu,good_cells);
%[a,b]=CCG(sp.st(idx),double(sp.clu(idx)),'binSize',[0.001],'duration',[0.5]);

%mex -O CCGHeart.c
spikes = double(sp.clu(idx));
tempSP=[ones(size(spikes)) spikes spikes];
mono=bz_MonoSynConvClick(double(tempSP),sp.st(idx),'plot',true);


%% VERIFY
    MM_animal=aggregateData.avgMM(aggregateData.AID==3 & aggregateData.CGS==2,:);
    MM_animal=MM_animal-mean(MM_animal(:,380:399),2);

for ii=1:length(mono.sig_con)
    figure
    subplot(2,1,1)
    plot([-.1:0.0004:.1],mono.ccgR(:,mono.sig_con(ii,1),mono.sig_con(ii,2)));
    grid on
    
    subplot(2,1,2)
    tvec=[-4:0.02:4];
    plot(tvec,MM_animal(good_cells==mono.sig_con(ii,1),:));
    hold on
    plot(tvec,MM_animal(good_cells==mono.sig_con(ii,2),:));
    xlim([-.5 2])
    grid on
    legend({'Pre','Post'})
    xlabel(['time'])
    ylabel(['firing rate change'])
    title('Mismatch Response')
end
%%
Depth_animal = aggregateData.DEPTH(aggregateData.AID==3 & aggregateData.CGS==2,:);
d_pre=[];
d_post=[];
for ii=1:length(mono.sig_con)
    d_pre(end+1)=Depth_animal(good_cells==mono.sig_con(ii,1));
    d_post(end+1) = Depth_animal(good_cells==mono.sig_con(ii,2));
    
end
figure
plot([d_pre',d_post']')
set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'})
xlim([.9 2.1])
ylabel('Depth')