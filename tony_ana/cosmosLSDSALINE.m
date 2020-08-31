%%
% salineData.GratAmpPreSaline=GratAmpPreSaline;
% salineData.GratAmpPostSaline=GratAmpPostSaline;
% salineData.cell_hemis = cell_hemis;
% salineData.cell_regions=cell_regions;
% salineData.GratSalinePre = GratSalinePre;
% salineData.GratSalinePost = GratSalinePost;
%% plot Saline Data

load('salineData')
frameRate = 1/0.034;

thresh = 10;
responsiveIDX = any(salineData.GratAmpPreSaline>thresh,2) & any(salineData.GratAmpPostSaline>thresh,2);

AllGratingsPre = mean(salineData.GratSalinePre,3);
AllGratingsPost = mean(salineData.GratSalinePost,3);

IDX = ismember(salineData.cell_regions,'VIS') & responsiveIDX' & logical(salineData.cell_hemis)';
figure
params = struct();
params.winIDX = (1:size(AllGratingsPre,2))-89; % in sample number, has to be the same as number of frames in resp
params.masterTime = params.winIDX/frameRate; % convert sample number to time
params.xLim = [-1 8]; % shows only from -1s to 3s
plotAVGSEM(AllGratingsPre(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89)
params.winIDX = (1:size(AllGratingsPost,2))-89; % in sample number, has to be the same as number of frames in resp
params.masterTime = params.winIDX/frameRate; % convert sample number to time
plotAVGSEM(AllGratingsPost(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89,'col',[1 0 0])
title('Saline')
legend({'Pre','Post'})
xlabel(nnz(IDX))

%%
load('lsdData')
frameRate = 1/0.034;

responsiveIDX = any(lsdData.GratAmpPrelsd>thresh,2) & any(lsdData.GratAmpPostlsd>thresh,2);

AllGratingsPre = mean(lsdData.GratlsdPre,3);
AllGratingsPost = mean(lsdData.GratlsdPost,3);

IDX = ismember(lsdData.cell_regions,'VIS') & responsiveIDX' & logical(lsdData.cell_hemis)';
figure
params.winIDX = (1:size(AllGratingsPre,2))-89; % in sample number, has to be the same as number of frames in resp
params.masterTime = params.winIDX/frameRate; % convert sample number to time
plotAVGSEM(AllGratingsPre(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89)
params.winIDX = (1:size(AllGratingsPost,2))-89; % in sample number, has to be the same as number of frames in resp
params.masterTime = params.winIDX/frameRate; % convert sample number to time
plotAVGSEM(AllGratingsPost(IDX,:)',gca,'parameters',params,'ms',true,'baseline',59:89,'col',[1 0 0])
title('LSD')
legend({'Pre','Post'})
xlabel(nnz(IDX))
%%
AllGratingsPreSaline = mean(salineData.GratSalinePre,3);
AllGratingsPostSaline = mean(salineData.GratSalinePost,3);
AGPSR = mean(AllGratingsPreSaline(:,95:155),2)-mean(AllGratingsPreSaline(:,39:89),2);

figure('Position',[440    59   460   739])
IDX = ismember(salineData.cell_regions,'VIS') & logical(salineData.cell_hemis)';
[~,sidx]=sort(AGPSR(IDX),'descend');
clim=[-20 20];
subplot(2,2,1)
tmp = bsxfun(@minus,AllGratingsPreSaline(IDX,:),mean(AllGratingsPreSaline(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('Saline Pre Hemi 1')
subplot(2,2,2)
tmp = bsxfun(@minus,AllGratingsPostSaline(IDX,:),mean(AllGratingsPostSaline(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('Saline Post Hemi 1')

IDX = ismember(salineData.cell_regions,'VIS') & ~logical(salineData.cell_hemis)';
[~,sidx]=sort(AGPSR(IDX),'descend');
clim=[-20 20];
subplot(2,2,3)
tmp = bsxfun(@minus,AllGratingsPreSaline(IDX,:),mean(AllGratingsPreSaline(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('Saline Pre Hemi 0')
subplot(2,2,4)
tmp = bsxfun(@minus,AllGratingsPostSaline(IDX,:),mean(AllGratingsPostSaline(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('Saline Post Hemi 0')
%%

AllGratingsPrelsd = mean(lsdData.GratlsdPre,3);
AllGratingsPostlsd = mean(lsdData.GratlsdPost,3);
AGPSR = mean(AllGratingsPrelsd(:,95:155),2)-mean(AllGratingsPrelsd(:,39:89),2);

figure('Position',[440    59   460   739])
IDX = ismember(lsdData.cell_regions,'VIS') & logical(lsdData.cell_hemis)';
[~,sidx]=sort(AGPSR(IDX),'descend');
clim=[-20 20];
subplot(2,2,1)
tmp = bsxfun(@minus,AllGratingsPrelsd(IDX,:),mean(AllGratingsPrelsd(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('lsd Pre Hemi 1')
subplot(2,2,2)
tmp = bsxfun(@minus,AllGratingsPostlsd(IDX,:),mean(AllGratingsPostlsd(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('lsd Post Hemi 1')

IDX = ismember(lsdData.cell_regions,'VIS') & ~logical(lsdData.cell_hemis)';
[~,sidx]=sort(AGPSR(IDX),'descend');
clim=[-20 20];
subplot(2,2,3)
tmp = bsxfun(@minus,AllGratingsPrelsd(IDX,:),mean(AllGratingsPrelsd(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('lsd Pre Hemi 0')
subplot(2,2,4)
tmp = bsxfun(@minus,AllGratingsPostlsd(IDX,:),mean(AllGratingsPostlsd(IDX,39:89),2));

imagesc(tmp(sidx,:),clim)
title('lsd Post Hemi 0')