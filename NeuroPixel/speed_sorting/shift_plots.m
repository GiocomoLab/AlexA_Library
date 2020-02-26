FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
DEPTH=[];
reg= [];
order = [1,2];
for iF=1:size(output,1)
    for iRep=1:2
    if isempty(output{iF}{iRep})
        disp('help')
        continue
    end
 
    FAST =cat(1,FAST,output{iF}{iRep}.all_fast);
    SLOW = cat(1,SLOW,output{iF}{iRep}.all_slow);
    GAIN = cat(1,GAIN,output{iF}{iRep}.all_gain);
    STAB = cat(1,STAB,output{iF}{iRep}.similarity);
    FACT = cat(1,FACT,output{iF}{order(iRep)}.factors');
    DEPTH = cat(1,DEPTH,output{iF}{iRep}.depth');
    reg = cat(2,reg,output{iF}{iRep}.region);
    end
end


figure
idx = STAB>.4 & ismember(reg,'VISp')';
med_fact = nanmedian(FACT(idx));

idx_1 = idx & FACT<-.1;
idx_2 = idx & FACT>-0.05;
subplot(2,1,1)
plot(x_vec,nanmean(FAST(idx_1,:)))
hold on
plot(x_vec,nanmean(SLOW(idx_1,:)));
plot(x_vec,nanmean(GAIN(idx_1,:)))
legend({'fast','slow','gain'})
grid on
box off
title('delayed cells')
subplot(2,1,2)
plot(x_vec,nanmean(FAST(idx_2,:)))
hold on
plot(x_vec,nanmean(SLOW(idx_2,:)));
plot(x_vec,nanmean(GAIN(idx_2,:)))
legend({'fast','slow','gain'})
grid on
box off
title('non delayed cells')
%% mean across all cells
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
DEPTH = [];
reg= [];
DCORR=[];
BLCORR = [];
for iF=1:size(output,1)
    
    
    for iRep = 1:2
        if ~isempty(output{iF}{iRep})
            FAST =cat(1,FAST,output{iF}{iRep}.all_fast);
            SLOW = cat(1,SLOW,output{iF}{iRep}.all_slow);
            GAIN = cat(1,GAIN,output{iF}{iRep}.all_gain);
            STAB = cat(1,STAB,output{iF}{iRep}.stability);
            FACT = cat(1,FACT,output{iF}{iRep}.factors');
            reg = cat(2,reg,output{iF}{iRep}.region);
            tmp = output{iF}{iRep}.correlation_shifted - output{iF}{iRep}.correlation_noshift;
            DCORR = cat(1,DCORR,(tmp'));
            BLCORR = cat(1,BLCORR,output{iF}{iRep}.correlation_noshift');
            tmp = output{iF}{iRep}.depth;
            if isrow(tmp)
                tmp = tmp';
            end
            DEPTH = cat(1,DEPTH,tmp);

        end
    end
end


figure
region = 'VISp';
idx = STAB>.5 & startsWith(reg,region)';
idx = DCORR>0 & startsWith(reg,region)' & STAB>.4;
boundedline(x_vec,nanmean(FAST(idx,:)),nanstd(FAST(idx,:))/sqrt(size(FAST(idx,:),1)),'alpha','cmap',[0 0 0])
hold on
%plot(x_vec,nanmean(FAST))
boundedline(x_vec,nanmean(SLOW(idx,:)),nanstd(SLOW(idx,:))/sqrt(size(FAST(idx,:),1)),'alpha')
boundedline(x_vec,nanmean(GAIN(idx,:)),nanstd(GAIN(idx,:))/sqrt(size(FAST(idx,:),1)),'alpha','cmap',get_color(gain,100))
legend({'fast','slow','gain'})
set(gcf,'Renderer','Painters')
grid on
box off
f=figure();
f.Name = region;
idx = STAB>0.4 & startsWith(reg,region)';
idx2 = STAB>0.4 & startsWith(reg,region)' & DCORR>0;

subplot(1,3,1)
scatter(DCORR(idx),DEPTH(idx),24,'.')
hold on
scatter(DCORR(idx2),DEPTH(idx2),15,[1 0 0 ],'o') 
set(gca,'YDir','reverse')
title('Diff in Corr')
ylabel('Depth')
subplot(1,3,2)
scatter(BLCORR(idx),DEPTH(idx),24,'.')
hold on
scatter(BLCORR(idx2),DEPTH(idx2),15,[1 0 0],'o')
set(gca,'YDir','reverse')
title('BL Corr')
ylabel('Depth')
subplot(1,3,3)
scatter(FACT(idx),DEPTH(idx),24,'.')
hold on
scatter(FACT(idx2),DEPTH(idx2),15,[1 0 0],'o')

set(gca,'YDir','reverse')
title('Delay Factor')
ylabel('Depth')
figure
subplot(1,2,1)
scatter(BLCORR(idx),DCORR(idx))
xlabel('BL Corr')
ylabel('Corr Diff')
subplot(1,2,2)
scatter(FACT(idx),DCORR(idx))
xlabel('Factor')
ylabel('Diff')
%saveas(gcf,sprintf('/oak/stanford/groups/giocomo/attialex/FIGURES/peaks_%s.pdf',region))
%% mean across site
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
DCORR = [];
reg= [];
region = 'VISp';

for iF=1:size(output,1)
    if isempty(output{iF})
        continue
    end
    if numel(output{iF}{1}.stability) ~= numel(output{iF}{1}.region)
        disp(iF)
        continue
    end
    for iRep =1:2
        if ~isempty(output{iF}{iRep})
            idx_region = startsWith(output{iF}{iRep}.region,region)';
            idx = output{iF}{iRep}.similarity>.4 & startsWith(output{iF}{iRep}.region,region)';
            if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.2
                sel_reg = mean(output{iF}{iRep}.similarity(idx));

                FAST =cat(1,FAST,nanmean(output{iF}{iRep}.all_fast(idx,:),1));
                SLOW = cat(1,SLOW,nanmean(output{iF}{iRep}.all_slow(idx,:),1));
                GAIN = cat(1,GAIN,nanmean(output{iF}{iRep}.all_gain(idx,:),1));
                tmp = output{iF}{iRep}.correlation_shifted - output{iF}{iRep}.correlation_noshift;
                DCORR = cat(1,DCORR,nanmean(tmp(idx)));
            end
        end
    end
end


figure('Position',[1356         405         475         525])
subplot(4,1,[1 2 3])
boundedline(x_vec,nanmean(FAST),nanstd(FAST)/sqrt(size(FAST,1)),'alpha','cmap',[0 0 0])
hold on
%plot(x_vec,nanmean(FAST))
boundedline(x_vec,nanmean(SLOW),nanstd(SLOW)/sqrt(size(FAST,1)),'alpha')
boundedline(x_vec,nanmean(GAIN),nanstd(GAIN)/sqrt(size(FAST,1)),'alpha','cmap',get_color(gain,100))
title(region)

legend({'fast','slow','gain'})
set(gcf,'Renderer','Painters')
grid on
box off
subplot(4,1,4)
plotSpread(DCORR,'xyOri','flipped')
title('Change in Correlation testset')
saveas(gcf,sprintf('/oak/stanford/groups/giocomo/attialex/FIGURES/peaks_%s_new.pdf',region))

%%
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
DCORR = [];
reg= [];
region = 'VISp6';

for iF=1:size(output,1)
    if isempty(output{iF})
        continue
    end
    if numel(output{iF}{1}.stability) ~= numel(output{iF}{1}.region)
        disp(iF)
        continue
    end
    for iRep =1:2
        if ~isempty(output{iF}{iRep})
            idx_region = startsWith(output{iF}{iRep}.subregion,region)';
            idx = output{iF}{iRep}.similarity>.4 & startsWith(output{iF}{iRep}.subregion,region)';
            if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.2
                sel_reg = mean(output{iF}{iRep}.similarity(idx));

                FAST =cat(1,FAST,nanmean(output{iF}{iRep}.all_fast(idx,:),1));
                SLOW = cat(1,SLOW,nanmean(output{iF}{iRep}.all_slow(idx,:),1));
                GAIN = cat(1,GAIN,nanmean(output{iF}{iRep}.all_gain(idx,:),1));
                tmp = output{iF}{iRep}.correlation_shifted - output{iF}{iRep}.correlation_noshift;
                DCORR = cat(1,DCORR,nanmean(tmp(idx)));
            end
        end
    end
end


figure('Position',[1356         405         475         525])
subplot(4,1,[1 2 3])
boundedline(x_vec,nanmean(FAST),nanstd(FAST)/sqrt(size(FAST,1)),'alpha','cmap',[0 0 0])
hold on
%plot(x_vec,nanmean(FAST))
boundedline(x_vec,nanmean(SLOW),nanstd(FAST)/sqrt(size(FAST,1)),'alpha')
boundedline(x_vec,nanmean(GAIN),nanstd(FAST)/sqrt(size(FAST,1)),'alpha','cmap',get_color(gain,100))
title(region)

legend({'fast','slow','gain'})
set(gcf,'Renderer','Painters')
grid on
box off
subplot(4,1,4)
plotSpread(DCORR,'xyOri','flipped')
title('Change in Correlation testset')
%saveas(gcf,sprintf('/oak/stanford/groups/giocomo/attialex/FIGURES/peaks_%s_new.pdf',region))
