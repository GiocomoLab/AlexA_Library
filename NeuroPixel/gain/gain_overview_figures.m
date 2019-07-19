%% normalized by c=100
SIM_ALL = [];
for iF=1:length(MERGED)
g_values = unique(MERGED(iF).gain);
sidx = MERGED(iF).stability;
[a,b]=sort(sidx,'descend');
frac = 1;
n=floor(numel(b)*frac);
tmp = [];
avgSimilarity=ones(n*numel(g_values),2);

for iC=numel(g_values):-1:1
    onsets = strfind(MERGED(iF).gain' == g_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity(:,b(1:n))',onsets,'win',[0 3]);
    avg = mean(snps,3);
    
    start = (iC-1)*n+1;
    stop = iC*n;
    if g_values(iC)==1
        tmp = mean(avg,2);
        avgSimilarity(start:stop,1)=mean(avg,2);
        avgSimilarity(start:stop,2)=g_values(iC)*ones(n,1);
    else
        
        avgSimilarity(start:stop,1)=mean(avg,2)./tmp;
        avgSimilarity(start:stop,2)=g_values(iC)*ones(n,1);
    end
end
% figure
% boxplot(avgSimilarity(:,1),avgSimilarity(:,2));
SIM_ALL=cat(1,SIM_ALL,avgSimilarity);
end
figure
 boxplot(SIM_ALL(:,1),SIM_ALL(:,2))
 xlim([0.5000    6.5])
 ylim([-2 2])
 title(sprintf('Normalized, Frac: %.2f',frac))
 xlabel('Gain')
 %% not normalized
 SIM_ALL = [];
for iF=1:length(MERGED)
g_values = unique(MERGED(iF).gain);
sidx = MERGED(iF).stability;
[a,b]=sort(sidx,'descend');
frac = 1;
n=floor(numel(b)*frac);
tmp = [];
avgSimilarity=ones(n*numel(g_values),2);

for iC=numel(g_values):-1:1
    onsets = strfind(MERGED(iF).gain' == g_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity(:,b(1:n))',onsets,'win',[0 3]);
    avg = mean(snps,3);
    
    start = (iC-1)*n+1;
    stop = iC*n;
    
    
        
        avgSimilarity(start:stop,1)=mean(avg,2);
        avgSimilarity(start:stop,2)=g_values(iC)*ones(n,1);

end
% figure
% boxplot(avgSimilarity(:,1),avgSimilarity(:,2));
SIM_ALL=cat(1,SIM_ALL,avgSimilarity);
end
figure
 boxplot(SIM_ALL(:,1),SIM_ALL(:,2))
 xlim([0.5000    7.5])
 ylim([-.5 1])
 title(sprintf('Frac: %.2f',frac))
 xlabel('gain')
 
  %% not normalized, averaged across animals
 g_values = [0.5 0.6 0.7 0.8];
  SIM_ALL = zeros(15,numel(g_values));

for iF=1:length(MERGED)
sidx = MERGED(iF).stability;
[a,b]=sort(sidx,'descend');
frac = .2;
n=floor(numel(b)*frac);
tmp = [];

for iC=numel(g_values):-1:1
    onsets = strfind(MERGED(iF).gain' == g_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity(:,b(1:n))',onsets,'win',[0 9]);
    avg = mean(snps,3);
    
    avgPerCell = mean(avg,2);
    avgSite=mean(avgPerCell);
    
        
    SIM_ALL(iF,iC)=avgSite;
      

end
% figure
% boxplot(avgSimilarity(:,1),avgSimilarity(:,2));
end
figure
plot(SIM_ALL','x')
hold on
plot(nanmean(SIM_ALL,1),'ok','MarkerSize',15)
set(gca,'XTick',1:numel(g_values),'XTickLabel',g_values)
 title(sprintf('Frac: %.2f',frac))
xlabel('Contrast')

SIM_norm = bsxfun(@rdivide,SIM_ALL,SIM_ALL(:,end));
figure
plot(SIM_norm','x')
hold on
plot(nanmean(SIM_norm,1),'ok','MarkerSize',15)
set(gca,'XTickLabel',g_values)
 title(sprintf('Frac: %.2f',frac))
xlabel('Contrast')
%% collapse maps to one map per 4 stimuli
for iF=2:15
bl_idx = MERGED(iF).contrast==100 & MERGED(iF).gain==1;
bl_map = nanmean(MERGED(iF).maps_smoothed(:,:,bl_idx),3);
stepsize=4;
n_steps = floor(numel(MERGED(iF).contrast)/stepsize);
SIM_ALL = zeros(numel(MERGED(iF).stability),n_steps);
G_ALL = zeros(1,n_steps);
for steps = 1:n_steps
    batch = (steps-1)*stepsize+1:steps*stepsize;
    if numel(unique(MERGED(iF).gain(batch)))>1
        error('multiple contrast values in this step')
    end
    batch_map = nanmean(MERGED(iF).maps_smoothed(:,:,batch),3);
    batch_sim = corr(batch_map',bl_map');
    batch_sim = diag(batch_sim);
    batch_gain = unique(MERGED(iF).gain(batch));
    SIM_ALL(:,steps)=batch_sim;
    G_ALL(steps)=batch_gain;
end
[~,sidx]=sort(MERGED(iF).stability,'descend');
figure
imagesc(SIM_ALL(sidx,:))
set(gca,'XTick',1:numel(G_ALL),'XTickLabel',G_ALL)
xlim([1.5 22.5])

end
%%
%% collapse maps to one map over 10 trials
%%
for iF=2:15
    figure
    imagesc(MERGED(iF).avgMapsSimilarity)
end

%%

for iF=1:5
    figure('Position',[24          42        1562         936])
    [~,sidx]=sort(MERGED(iF).stability,'descend');
    for iC=1:numel(MERGED(iF).avgMaps)
        subplot(1,numel(MERGED(iF).avgMaps),iC)
        normMap = bsxfun(@rdivide,MERGED(iF).avgMaps{iC}(sidx(1:20),:),sum(MERGED(iF).avgMaps{end}(sidx(1:20),:),2));
        imagesc(normMap)
        set(gca,'YTick',[])
    end
end
%%
col = winter(7);
for iF=2:15
        [~,sidx]=sort(MERGED(iF).stability,'descend');

figure;
    ha=tight_subplot(4,4,.01);


for iCell = 1:16
    axes(ha(iCell))
for iC=1:numel(MERGED(iF).avgMaps)
    hold on
    
    if iC==1
        plot(MERGED(iF).avgMaps{iC}(sidx(iCell),:)*50,'Color',[0 0 0])
    else
        plot(MERGED(iF).avgMaps{iC}(sidx(iCell),:)*50,'Color',col(iC,:))
    end
end
end
end
legend(array2stringCell(unique(MERGED(iF).contrast)))



%%

g_values = unique(MERGED(1).gain);
snps_win=[-4 7];
ALL_SIM = zeros(numel(g_values),numel(snps_win(1):snps_win(2)),numel(MERGED));

for iF=1:numel(MERGED)
    sidx = MERGED(iF).stability;
[~,sidx]=sort(sidx,'descend');
frac = .5;
n=floor(numel(sidx)*frac);
avgSimilarity=ones(numel(sidx)*numel(g_values),2);
for iC=1:numel(g_values)
    onsets = strfind(MERGED(iF).gain' == g_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity',onsets,'win',snps_win);
    avg = mean(snps,3);
    avg_all= mean(avg);
    avg_stab = mean(avg(sidx(1:n),:));
    
    
    ALL_SIM(iC,:,iF)=avg_stab;
    %subplot(1,3,3)
    
    
    
end
end
figure
plot(mean(ALL_SIM,3)')

figure
for iC=1:numel(g_values)
subplot(1,numel(g_values),iC)
plot(squeeze(ALL_SIM(iC,:,:)));
end

figure
col = winter(numel(g_values));
hold on
for iC=1:numel(g_values)-1
errorbar(snps_win(1):snps_win(2),nanmean(ALL_SIM(iC,:,:),3)',nanstd(ALL_SIM(iC,:,:),[],3)'/sqrt(numel(MERGED)),'Color',col(iC,:))
end
legend(array2stringCell(g_values'),'Location','southeast')
patch([ -.1 3.1  3.1 -.1],[0 0 1 1],[245/255, 224/255, 66/255],'FaceAlpha',.2,'EdgeColor','None')
ylim([0 .7])
title(sprintf('Fraction of neurons: %.2f',frac))
% 
% figure
% plot(squeeze(ALL_SIM(7,:,:)));

