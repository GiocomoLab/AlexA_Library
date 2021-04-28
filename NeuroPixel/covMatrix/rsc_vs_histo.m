matfiles = dir('/Volumes/T7/attialex/tbtxcorr_rsc_0.5/*.mat');
table_path = '/Users/attialex/code/AlexA_Library/NeuroPixel/histology/rsc_tables';
SHIFT = [];
STAB = [];
COORDS = [];
REGION = {};
DEPTH = [];
for iF=1:numel(matfiles)
    data = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    idx = startsWith(data.region,'RS');
    clus = data.good_cells(idx);
    SHIFT = cat(1,SHIFT,data.shiftMat(idx,:,:));
    STAB = cat(1,STAB,data.corrMat(idx,:,:));
    REGION = cat(2,REGION,data.region(idx));
    if isfield(data.anatomy,'depth_shifted')
        tmp_depth = data.anatomy.depth_shifted;
    else
        tmp_depth = data.anatomy.depth;
    end
    
    DEPTH = cat(2,DEPTH,tmp_depth(idx));
    sn = strrep(matfiles(iF).name,'_1.mat','_anatomy.csv');
    sn = [matfiles(iF).name(1:end-6),'_anatomy.csv'];
    data_table = readtable(fullfile(table_path,sn));
    nClu = numel(clus);
    tmp_coords = zeros(nClu,3);
    for iClu =1:nClu
        clu_idx = find(data_table.cluster_id==clus(iClu));
        tmp_coords(iClu,1) = data_table.xpos(clu_idx);
        tmp_coords(iClu,2) = data_table.ypos(clu_idx);
        tmp_coords(iClu,3) = data_table.zpos(clu_idx);
    end
    COORDS = cat(1,COORDS,tmp_coords);
end
%%
ss = nanmean(nanmean(STAB(:,1:6,1:6),2),3);

figure
tmp_d = DEPTH;
%tmp_d = COORDS(:,3)/1000;
tmp_d(ss<=0.5)=nan;
edges = prctile(tmp_d,[0 33 66 100]);
tmp = discretize(tmp_d,edges);
for ii=1:numel(edges)-1
    subplot(numel(edges)-1,2,(ii-1)*2+1)
    hold on
    idx = tmp==ii;
    imagesc(1:16,1:16,squeeze(nanmean(STAB(idx,:,:))),[0,0.7])
    title(sprintf('n=%d, d=%.3f',nnz(idx),nanmean(tmp_d(idx))))
    axis image
    subplot(numel(edges)-1,2,(ii)*2)
    hold on
    imagesc(1:16,1:16,squeeze(nanmean(SHIFT(idx,:,:))),[-5 5])
    axis image
end
%%
%% split by region and AP
[rr,~,bb]=unique(REGION);
X = COORDS(:,1);
figure
hold on
cols = lines(3);
for ii=1:3
    
    idx = ss>0.5 & bb==ii;
    subplot(1,2,1)
    hold on
remap = nanmean(nanmean(STAB(idx,1:6,7:10),2),3);
scatter(X(idx),remap,45,cols(ii,:),'.')

    subplot(1,2,2)
    hold on
shift = nanmean(nanmean(SHIFT(idx,1:6,7:10),2),3);
scatter(X(idx),shift,45,cols(ii,:),'.')
end
subplot(1,2,1)
legend(rr)
xlabel('A -> P')
ylabel('stability relative to baseline')

subplot(1,2,2)
legend(rr)
xlabel('A -> P')
ylabel('map shift')

%% split by region and depth
[rr,~,bb]=unique(REGION);

figure
hold on
cols = lines(3);
for ii=1:3
    
    idx = ss>0.5 & bb==ii;
    subplot(1,2,1)
    hold on
remap = nanmean(nanmean(STAB(idx,1:6,7:10),2),3);
scatter(DEPTH(idx),remap,45,cols(ii,:),'.')

    subplot(1,2,2)
    hold on
shift = nanmean(nanmean(SHIFT(idx,1:6,7:10),2),3);
scatter(DEPTH(idx),shift,45,cols(ii,:),'.')
end
subplot(1,2,1)
legend(rr)
xlabel('distance from surface')
ylabel('stability relative to baseline')

subplot(1,2,2)
legend(rr)
xlabel('distance from surface')
ylabel('map shift')
%%

ss = nanmean(nanmean(STAB(:,1:6,1:6),2),3);

figure
[rr,~,bb]=unique(REGION);
for ii=1:length(rr)
    subplot(numel(rr),2,(ii-1)*2+1)
    hold on
    idx = bb==ii & ss>0.5;
    imagesc(1:16,1:16,squeeze(nanmean(STAB(idx,:,:))),[0,0.7])
    title(sprintf('n=%d',nnz(idx)))
    axis image
    subplot(numel(rr),2,(ii)*2)
    hold on
    imagesc(1:16,1:16,squeeze(nanmean(SHIFT(idx,:,:))),[-5 5])
    title(rr{ii})
    axis image
end
    