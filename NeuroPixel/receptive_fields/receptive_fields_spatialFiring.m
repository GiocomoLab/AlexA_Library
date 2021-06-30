data_table = readtable('/Volumes/T7/attialex/receptive_fields/vi_trippy_pairs.xlsx');

trippy_path = 'F:\Alex\receptive_fields';
trippy_path = '/volumes/T7/attialex/receptive_fields';
data_path = '/volumes/T7/attialex/NP_Data_corrected';
opt = load_default_opt;
%% grouped by recording on same animal and hemisphere
groups=unique(data_table.group);
marker = {'x','.','o'};
cols = brewermap(3,'Set1');
for iG=1%:numel(groups)
    idx = find(data_table.group == groups(iG));
    animal = data_table.Animal(idx(1));
    probe_nbrs = data_table.Probe_nbr(idx);
    if length(idx)<2
        continue
    end
    
    TC=cell(1,numel(idx));
    RF = cell(1,numel(idx));
    MU = cell(1,numel(idx));
    for ii=1:numel(idx)
        data = load(fullfile(trippy_path,data_table.Trippy_name{idx(ii)},'receptive_fields.mat'));
        
        rf = [];
        rf_noT = [];
        cells_with_rf=[];
        allMU=[];
        for iC=1:numel(data.fields)
            p=data.fields{iC};
            if length(p)>0
                ff= squeeze(data.staMat(:,:,3:7,iC));
                tmp = mean(ff,3);
                tmp = imgaussfilt(tmp);
                Z = zscore(tmp,[],'All');
                rf_noT = cat(3,rf_noT,Z);
                Z(abs(Z)<2.5)=0;
                rf = cat(3,rf,Z);
                cells_with_rf(end+1)=data.good_cells(iC);
                tmpMU=[];
                for iR=1:numel(p)
                    tmpMU=cat(1,tmpMU,p{iR}.mu);
                end
                if size(tmpMU,1)>1
                    tmpMU = mean(tmpMU);
                end
                
                allMU=cat(1,allMU,tmpMU);
                
                
                
            end
            
        end
        spatial_data = load(fullfile(data_path,data_table.Spatial_name{idx(ii)}));
        [corrMat,frMat,~]=trialCorrMat(cells_with_rf,5:20,spatial_data,opt);
        tc = squeeze(mean(frMat,2));
        tcZ = zscore(tc,[],2);
        pdist_tc = pdist(tcZ,'Correlation');
        
        pdist_rf = pdist(reshape(rf_noT,[],size(rf,3))','Correlation');
        v_idx = isfinite(pdist_rf) & isfinite(pdist_tc);
        [a,b]=corrcoef(pdist_rf(v_idx),pdist_tc(v_idx))
        TC{ii}=tc;
        RF{ii}=rf_noT;
        if startsWith(data_table.Hemisphere{idx(ii)},'L')
            allMU(:,1)=176-allMU(:,1);
        end
        MU{ii}=allMU;
    end
    figure('Color','White')
    subplot(1,3,1)
    hold on
    for ii=1:numel(idx)
        scatter(MU{ii}(:,1),MU{ii}(:,2),45,cols(ii,:),'.')
    end
    set(gca,'YDir','reverse')
xlim([0,176])
ylim([0,96])
xlabel('L -> M')
ylabel('D ->U')
title('Individual Receptive Field Centers')
    allMU = cat(1,MU{:});
    allTC = cat(1,TC{:});
    [Y,E]=discretize(allMU(:,1),3);
    
    subplot(1,3,2)
    [~,sid]=sort(allMU(:,1));
    imagesc(zscore(allTC(sid,:),[],2),[-2 2])
    title(groups(iG))
    
    subplot(1,3,3)
    [~,sid]=sort(allMU(:,2));
    imagesc(zscore(allTC(sid,:),[],2),[-2 2])
    title(sprintf('Hemisphere %s',data_table.Hemisphere{idx(1)}))
end
%%
figure('Render','Painters','Color','White','Position',[441   575   735   224])
[~,sid]=sort(allMU(:,1));
    imagesc(1:2:400,1:numel(sid),zscore(allTC(sid,:),[],2),[-2 2])
    colormap(brewermap(20,'BuPu'))
xlabel('Position [cm]')
saveas(gcf,'/users/attialex/desktop/spatial_maps_group_1.pdf')
%%
[a,sid]=sort(pdist_tc);
[ii,jj]=triind2sub([size(rf,3),size(rf,3)],sid);
figure
for iR=1:15
    subplot(1,3,1)
    imagesc(squeeze(rf_noT(:,:,ii(iR))),[-5 5])
    subplot(1,3,2)
    imagesc(squeeze(rf_noT(:,:,jj(iR))),[-5 5])
    subplot(1,3,3)
    plot(tc(ii(iR),:))
    hold on
    plot(tc(jj(iR),:))
    pause
    
    clf
    
    
end
%%
[a,sid]=sort(pdist_tc);
[ii,jj]=triind2sub([size(rf,3),size(rf,3)],sid);
figure
for iR=1:15
    subplot(1,2,1)
    imagesc(squeeze(rf_noT(:,:,ii(iR))),[-5 5])
    subplot(1,2,2)
    imagesc(squeeze(rf(:,:,jj(iR))),[-5 5])
    pause
    clf
end
%%

groups=unique(data_table.Hemisphere);
marker = {'x','.','o'};
cols = brewermap(3,'Set1');
for iG=2%:numel(groups)
    idx = find(startsWith(data_table.Hemisphere, groups{iG}));
    
    
    TC=cell(1,numel(idx));
    RF = cell(1,numel(idx));
    MU = cell(1,numel(idx));
    for ii=1:numel(idx)
        data = load(fullfile(trippy_path,data_table.Trippy_name{idx(ii)},'receptive_fields.mat'));
        if ~isfile(fullfile(data_path,[data_table.Spatial_name{idx(ii)},'.mat']))
            disp('skip')
            continue
        end
        spatial_data = load(fullfile(data_path,data_table.Spatial_name{idx(ii)}));
        
        rf = [];
        rf_noT = [];
        cells_with_rf=[];
        allMU=[];
        for iC=1:numel(data.fields)
            p=data.fields{iC};
            if length(p)>0
                ff= squeeze(data.staMat(:,:,3:7,iC));
                tmp = mean(ff,3);
                tmp = imgaussfilt(tmp);
                Z = zscore(tmp,[],'All');
                rf_noT = cat(3,rf_noT,Z);
                Z(abs(Z)<2.5)=0;
                rf = cat(3,rf,Z);
                cells_with_rf(end+1)=data.good_cells(iC);
                tmpMU=[];
                for iR=1:numel(p)
                    tmpMU=cat(1,tmpMU,p{iR}.mu);
                end
                if size(tmpMU,1)>1
                    tmpMU = mean(tmpMU);
                end
                
                allMU=cat(1,allMU,tmpMU);
                
                
                
            end
            
        end
        [corrMat,frMat,~]=trialCorrMat(cells_with_rf,5:20,spatial_data,opt);
        tc = squeeze(mean(frMat,2));
        tcZ = zscore(tc,[],2);
        
        TC{ii}=tc;
        RF{ii}=rf_noT;
        MU{ii}=allMU;
    end
end
%%
allMU = cat(1,MU{:});
allTC = cat(1,TC{:});
figure
subplot(1,2,1)
[~,sid]=sort(allMU(:,1));
imagesc(zscore(allTC(sid,:),[],2),[-2 2])
subplot(1,2,2)
[~,sid]=sort(allMU(:,2));
imagesc(zscore(allTC(sid,:),[],2),[-2 2])
%% all RF's projected onto right side

marker = {'x','.','o'};
cols = brewermap(3,'Set1');
nS=size(data_table.Trippy_name,1);
TC=cell(1,nS);
    RF = cell(1,nS);
    MU = cell(1,nS);
    GC = cell(1,nS);
    STAB = cell(1,nS);
for iS=1:nS
    
    
    
    
    data = load(fullfile(trippy_path,data_table.Trippy_name{iS},'receptive_fields.mat'));
    if ~isfile(fullfile(data_path,[data_table.Spatial_name{iS},'.mat']))
        disp('skip')
        continue
    end
    spatial_data = load(fullfile(data_path,data_table.Spatial_name{iS}));
    
    rf = [];
    rf_noT = [];
    cells_with_rf=[];
    allMU=[];
    for iC=1:numel(data.fields)
        p=data.fields{iC};
        if length(p)>0
            ff= squeeze(data.staMat(:,:,3:7,iC));
            tmp = mean(ff,3);
            tmp = imgaussfilt(tmp);
            Z = zscore(tmp,[],'All');
            rf_noT = cat(3,rf_noT,Z);
            Z(abs(Z)<2.5)=0;
            rf = cat(3,rf,Z);
            cells_with_rf(end+1)=data.good_cells(iC);
            tmpMU=[];
            for iR=1:numel(p)
                tmpMU=cat(1,tmpMU,p{iR}.mu);
            end
            if size(tmpMU,1)>1
                tmpMU = mean(tmpMU);
            end
            
            allMU=cat(1,allMU,tmpMU);
            
        end
        
    end
    [corrMat,frMat,~]=trialCorrMat(cells_with_rf,5:20,spatial_data,opt);
    tc = squeeze(mean(frMat,2));
    tcZ = zscore(tc,[],2);
    STAB{iS}=nanmean(nanmean(corrMat,2),3);
    TC{iS}=tc;
    RF{iS}=rf_noT;
    if startsWith(data_table.Hemisphere{iS},'L')
        allMU(:,1)=172-allMU(:,1);
    end
    MU{iS}=allMU;
    GC{iS}=cells_with_rf;
    
end
%%
trippy_name = {};
spatial_name={};
mouse_name = {};
cluID = [];
allMU = [];
for iS=1:14
    if ~isempty(GC{iS})
        nC=numel(GC{iS});
    trippy_name = cat(1,trippy_name,repmat(data_table.Trippy_name(iS),nC,1));
    spatial_name = cat(1,spatial_name,repmat(data_table.Spatial_name(iS),nC,1));
    mouse_name = cat(1,mouse_name,repmat(data_table.Animal(iS),nC,1));
    cluID = cat(1,cluID,GC{iS}');
    allMU = cat(1,allMU,MU{iS});
    end
end
rec_x=allMU(:,1);
rec_y=allMU(:,2);
cell_table = table(trippy_name,spatial_name,mouse_name,cluID,rec_x,rec_y);
save('/Users/attialex/code/AlexA_Library/NeuroPixel/receptive_fields/receptive_field_table.mat','cell_table')
    
%%
figure('Color','White')
allMU = cat(1,MU{:});
allTC = cat(1,TC{:});

subplot(1,2,1)
[~,sid]=sort(allMU(:,1));
imagesc(zscore(allTC(sid,:),[],2),[-2 2])
xlabel('Position [cm]')
ylabel('Units sorted by ML')
subplot(1,2,2)
[~,sid]=sort(allMU(:,2),'descend');
imagesc(zscore(allTC(sid,:),[],2),[-2 2])
xlabel('Position [cm]')
ylabel('Units sorted by DV')
%%
stab_bl = cat(1,STAB{:});
figure('Color','White')
scatter(allMU(:,1),allMU(:,2),15,cat(1,stab_bl),'filled')
set(gca,'YDir','reverse','CLim',[0 0.5])
xlim([0,176])
ylim([0,96])
xlabel('L -> M')
ylabel('D ->U')
c = colorbar;
c.Label.String = 'BL Spatial Stability';

%%
edges = prctile(allMU(:,1),[0 25, 50,75,100]);
edges = prctile(allMU(:,1),[0 33, 66,100]);

d =discretize(allMU(:,1),edges);
figure('Color','white','Render','Painters')
hold on
%cols=brewermap(numel(edges)-1,'*RdYlBu')
cols = brewermap(numel(edges)-1,'Set1')
leg={};
for ie=1:numel(edges)-1
    idx = d==ie;
    mu = nanmean(zscore(allTC(idx,:),[],2));
    st = nanstd(zscore(allTC(idx,:),[],2))/sqrt(nnz(idx));
    boundedline(1:2:400,mu,st,'cmap',cols(ie,:),'alpha')
    %plot(1:2:400,nanmean(zscore(allTC(idx,:),[],2)),'Color',cols(ie,:),'LineWidth',2)
    leg{end+1}=sprintf('Bin %d/%d',ie,numel(edges)-1);
end
legend(leg)
xlabel('Position [cm]')
ylabel('Normalized FR [zscored]')