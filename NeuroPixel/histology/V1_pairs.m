[corrMat,frMat] = getSpatialMapsForRegion(data_1,'VISp',5:20);
corr_t=0.5;
stab = nanmean(nanmean(corrMat,2),3);
tc_1 = squeeze(nanmean(frMat(stab>corr_t,:,:),2));

[corrMat,frMat] = getSpatialMapsForRegion(data_2,'VISp',5:20);
corr_t=0.5;
stab = nanmean(nanmean(corrMat,2),3);
tc_2 = squeeze(nanmean(frMat(stab>corr_t,:,:),2));

%%
cc=corr(tc_1',tc_2');
[mv,ii]=max(cc,[],'all','linear');
[r,c]=ind2sub(size(cc),ii)
figure
plot(tc_1(r,:));
hold on
plot(tc_2(c,:));
%%

X = cat(1,tc_1,tc_2);

dist_mat = squareform(pdist(X,'cosine'));
dist_across = dist_mat(size(tc_1,1)+1:size(X,1),1:size(tc_1,1));

[mv,ii]=min(dist_across',[],'all','linear');
%[mv,ii]=max((1-dist_across'),[],'all','linear');
[r,c]=ind2sub(size(dist_across'),ii);
figure
plot(tc_1(r,:));
hold on
plot(tc_2(c,:));

%%
[~,sid]=sort(reshape(dist_across',1,[]),'ascend','MissingPlacement','last');
figure
for ii=1:21
    [r,c]=ind2sub(size(dist_across'),sid(ii));
    
    plot(tc_1(r,:));
hold on
plot(tc_2(c,:));
title(finddelay(tc_1(r,:),tc_2(c,:)))
pause
clf
end

%% pairs of V1
v1_data= { {'AA1_190726_gain_1.mat','AA1_190728_gaincontrast10_1.mat'};
    {'AA1_190729_gain_1.mat','AA1_190730_gaincontrast10_1'};
    {'AA2_190809_gaincontrast10_2.mat','AA2_190810_contrast_1.mat'};
    {'AA3_190801_gain_3.mat','AA3_190802_contrast_1'};
    {'AA3_190804_gain_1.mat','AA3_190805_gaincontrast10_1.mat'};
    {'AA4_190801_gain_1','AA4_190803_contrast_1'};
    {'AA4_190805_gaincontrast10_2.mat','AA4_190806_gain_1.mat'};
    {'AA50_190930_gaincontrast10_1.mat','AA50_191001_gain_1.mat'};
    }
