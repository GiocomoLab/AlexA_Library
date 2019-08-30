%% normalized by c=100
SIM_ALL = [];
for iF=1:length(MERGED)
c_values = unique(MERGED(iF).contrast);
sidx = MERGED(iF).stability;
[a,b]=sort(sidx,'descend');
frac = .3;
n=floor(numel(b)*frac);
tmp = [];
avgSimilarity=ones(n*numel(c_values),2);

for iC=numel(c_values):-1:1
    onsets = strfind(MERGED(iF).contrast' == c_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity(:,b(1:n))',onsets,'win',[0 9]);
    avg = mean(snps,3);
    
    start = (iC-1)*n+1;
    stop = iC*n;
    if c_values(iC)==100
        tmp = mean(avg,2);
        avgSimilarity(start:stop,1)=mean(avg,2);
        avgSimilarity(start:stop,2)=c_values(iC)*ones(n,1);
    else
        
        avgSimilarity(start:stop,1)=mean(avg,2)./tmp;
        avgSimilarity(start:stop,2)=c_values(iC)*ones(n,1);
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
 xlabel('Contrast')
 %% not normalized
 SIM_ALL = [];
for iF=1:length(MERGED)
c_values = unique(MERGED(iF).contrast);
sidx = MERGED(iF).stability;
[a,b]=sort(sidx,'descend');
frac = .3;
n=floor(numel(b)*frac);
tmp = [];
avgSimilarity=ones(n*numel(c_values),2);

for iC=numel(c_values):-1:1
    onsets = strfind(MERGED(iF).contrast' == c_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity(:,b(1:n))',onsets,'win',[0 9]);
    avg = mean(snps,3);
    
    start = (iC-1)*n+1;
    stop = iC*n;
    
        
        avgSimilarity(start:stop,1)=mean(avg,2);
        avgSimilarity(start:stop,2)=c_values(iC)*ones(n,1);

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
 xlabel('Contrast')
 
  %% not normalized, averaged across animals
 SIM_ALL = zeros(15,7);
 c_values = [0 2 5 10 20 50 100];
for iF=1:length(MERGED)
sidx = MERGED(iF).stability;
[a,b]=sort(sidx,'descend');
frac = .5;
n=floor(numel(b)*frac);
tmp = [];

for iC=numel(c_values):-1:1
    onsets = strfind(MERGED(iF).contrast' == c_values(iC),[0 1])+1;
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
plot(nanmean(SIM_ALL,1),'xk','MarkerSize',15)
set(gca,'XTickLabel',c_values)
 title(sprintf('Frac: %.2f',frac))
xlabel('Contrast')

SIM_norm = bsxfun(@rdivide,SIM_ALL,SIM_ALL(:,end));
figure
plot(SIM_norm(:,7:-1:1)','x')
hold on
plot(nanmean(SIM_norm(:,7:-1:1),1),'.k','MarkerSize',15)
set(gca,'XTickLabel',c_values(end:-1:1))
 title(sprintf('Frac: %.2f',frac))
xlabel('Contrast')
%% collapse maps to one map per 10 stimuli
for iF=2:15
bl_idx = MERGED(iF).contrast==100 & MERGED(iF).gain==1;
bl_map = nanmean(MERGED(iF).maps_smoothed(:,:,bl_idx),3);
stepsize=10;
n_steps = floor(numel(MERGED(iF).contrast)/stepsize);
SIM_ALL = zeros(numel(MERGED(iF).stability),n_steps);
C_ALL = zeros(1,n_steps);
for steps = 1:n_steps
    batch = (steps-1)*stepsize+1:steps*stepsize;
    if numel(unique(MERGED(iF).contrast(batch)))>1
        error('multiple contrast values in this step')
    end
    batch_map = nanmean(MERGED(iF).maps_smoothed(:,:,batch),3);
    batch_sim = corr(batch_map',bl_map');
    batch_sim = diag(batch_sim);
    batch_contrast = unique(MERGED(iF).contrast(batch));
    SIM_ALL(:,steps)=batch_sim;
    C_ALL(steps)=batch_contrast;
end
[~,sidx]=sort(MERGED(iF).stability,'descend');
figure
imagesc(SIM_ALL(sidx,:))
set(gca,'XTick',2:2:30,'XTickLabel',C_ALL(2:2:end))
xlim([1.5 22.5])

end
%%
%% collapse maps to one map over 10 trials
for iF=2:15
bl_idx = MERGED(iF).contrast==100 & MERGED(iF).gain==1;
bl_map = nanmean(MERGED(iF).maps_smoothed(:,:,bl_idx),3);
stepsize=10;
n_steps = floor(numel(MERGED(iF).contrast)/stepsize);
SIM_ALL = zeros(numel(MERGED(iF).stability),n_steps);
C_ALL = zeros(1,n_steps);
for steps = 1:n_steps
    batch = (steps-1)*stepsize+1:steps*stepsize;
    if numel(unique(MERGED(iF).contrast(batch)))>1
        error('multiple contrast values in this step')
    end
    batch_map = nanmean(MERGED(iF).maps_smoothed(:,:,batch),3);
    batch_sim = corr(batch_map',bl_map');
    batch_sim = diag(batch_sim);
    batch_contrast = unique(MERGED(iF).contrast(batch));
    SIM_ALL(:,steps)=batch_sim;
    C_ALL(steps)=batch_contrast;
end
[~,sidx]=sort(MERGED(iF).stability,'descend');
figure
imagesc(SIM_ALL(sidx,:))
end
%%
for iF=2:15
    figure
    imagesc(MERGED(iF).avgMapsSimilarity)
end

%%

for iF=7
    figure('Position',[24          42        1562         936])
    [~,sidx]=sort(MERGED(iF).stability,'descend');
    for iC=1:7
        subplot(1,7,iC)
        normMap = bsxfun(@rdivide,MERGED(iF).avgMaps{iC}(sidx(1:20),:),sum(MERGED(iF).avgMaps{7}(sidx(1:20),:),2));
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
for iC=1:7
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
%% one of the single cell examples
iF=13;
%fn = fullfile(root,Files(iF).name);
%contrast = load(fn);
        [~,sidx]=sort(MERGED(iF).stability,'descend');

good_cells = contrast.sp.cids(contrast.sp.cgs==2);
raster_map(contrast,good_cells(sidx(8)));
axis tight
ylabel('trial')
xlabel('Position')
%%
iF=11;
fn = fullfile(root,Files(iF).name);
contrast = load(fn);
        [~,sidx]=sort(MERGED(iF).stability,'descend');

good_cells = contrast.sp.cids(contrast.sp.cgs==2);
raster_map(contrast,good_cells(sidx(1)));
axis tight
ylabel('trial')

xlabel('Position')
%%
iF=7;
fn = fullfile(root,Files(iF).name);
contrast = load(fn);
        [~,sidx]=sort(MERGED(iF).stability,'descend');

good_cells = contrast.sp.cids(contrast.sp.cgs==2);
raster_map(contrast,good_cells(sidx(1)));
axis tight
ylabel('trial')
xlabel('Position')
%%

ALL_SIM = zeros(7,10,14);

for iF=2:15
    sidx = MERGED(iF).stability;
[~,sidx]=sort(sidx,'descend');
frac = .1;
n=floor(numel(sidx)*frac);
c_values = unique(MERGED(iF).contrast);
avgSimilarity=ones(numel(sidx)*numel(c_values),2);
for iC=1:numel(c_values)
    onsets = strfind(MERGED(iF).contrast' == c_values(iC),[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity',onsets,'win',[0 9]);
    avg = mean(snps,3);
    avg_all= mean(avg);
    avg_stab = mean(avg(sidx(1:n),:));
    
    
    ALL_SIM(iC,:,iF-1)=avg_stab;
    %subplot(1,3,3)
    
    
    
end
end
figure
plot(mean(ALL_SIM,3)')

figure
for iC=1:7
subplot(1,7,iC)
plot(squeeze(ALL_SIM(iC,:,:)));
end

figure
col = winter(7);
hold on
for iC=1:7
errorbar(1:10,mean(ALL_SIM(iC,:,:),3)',std(ALL_SIM(iC,:,:),[],3)'/sqrt(14),'Color',col(iC,:))
end

figure
plot(squeeze(ALL_SIM(7,:,:)));

%%
win = [-10 9];
ALL_SIM = zeros(1,sum(abs(win))+1,14);

for iF=2:15
    sidx = MERGED(iF).stability;
[~,sidx]=sort(sidx,'descend');
frac = .3;
n=floor(numel(sidx)*frac);
c_values = unique(MERGED(iF).contrast);
avgSimilarity=ones(numel(sidx)*numel(c_values),2);
for iC=1
    onsets = strfind(MERGED(iF).gain' == 0.5,[0 1])+1;
    snps = extract_snps(MERGED(iF).similarity',onsets,'win',win);
    avg = mean(snps,3);
    avg_all= mean(avg);
    avg_stab = mean(avg(sidx(1:n),:));
    
    
    ALL_SIM(iC,:,iF-1)=avg_stab;
    %subplot(1,3,3)
    
    
    
end
end


figure
col = winter(7);
hold on
for iC=1
errorbar(win(1):win(2),nanmean(ALL_SIM(iC,:,:),3)',nanstd(ALL_SIM(iC,:,:),[],3)'/sqrt(14),'Color',col(iC,:))
end
ylabel('similarity to baseline')
xlabel('trial')

figure
plot(win(1):win(2),squeeze(ALL_SIM(1,:,:)));
hold on
plot(win(1):win(2),nanmean(ALL_SIM,3),'k','LineWidth',2)
ylabel('similarity to baseline')
xlabel('trial')