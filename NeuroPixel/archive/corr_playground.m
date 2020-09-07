

%%
dist=[];
for kk=1:20
    [IDX,~,sumd]=kmeans(squeeze(spatialMap(1,:,:))',kk);
    dist(end+1)=sum(sumd);
end

figure;plot(dist)
[IDX]=kmeans(squeeze(spatialMap(1,:,:))',4);
[a,sortK]=sort(IDX);
figure
imagesc(squeeze(correlation_All(sortK,sortK,1)),[0 0.4])
%%
ffSorted = zeros(size(correlation_All));
for iC=1:size(spatialMap,1)
    [IDX]=kmeans(squeeze(spatialMap(iC,2:38,:))',4);
    [a,sortK]=sort(IDX);
    ffSorted(:,:,iC)=correlation_All(sortK,sortK,iC);
end
view_stack(cat(2,correlation_All,ffSorted))

%%

trial_color = zeros(1,max(trial),3);
baseline = [0,0,0];
contrast = [.4 .4 .4];
gain_1 = [0 1 1];
gain_2 = [1 0 1];
for ii = 1:size(trial_color,2)
    if exist('trial_contrast','var')
        if trial_contrast(ii)<100
            trial_color(1,ii,:)=contrast;
        end
    end
    if exist('trial_gain','var')    
        if trial_gain(ii)<1
            trial_color(1,ii,:)=gain_1;
        end
        if trial_gain(ii)<.7
            trial_color(1,ii,:)=gain_2;
        end
    end
end
trial_color=trial_color(1,sortK,:);
tc=permute(trial_color,[2 1 3]);
figure
subplot('Position',[0.1 0.1 .8 .8 ])
img=nanmean(correlation_All,3);
imagesc(img(sortK,sortK),[-.1 0.4])
    %xlabel(Files(iF).name,'Interpreter','None')

subplot('Position',[0.1 0.9 .8 0.05 ])
imagesc(trial_color)
set(gca,'XTick',[],'YTick',[])
subplot('Position',[0.05 0.1 .05 0.8 ])
imagesc(tc)
set(gca,'XTick',[],'YTick',[])
