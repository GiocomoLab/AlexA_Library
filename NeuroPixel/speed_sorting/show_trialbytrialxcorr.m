
%%
files = dir('F:\temp\xcorrSpeed\*.mat');

ALL_SLOPES = struct();
ALL_SESSIONS=struct();
ALL_CLUSTERS = struct();
for iF=1:numel(files)
    d=load(fullfile(files(iF).folder,files(iF).name));
    data_out = d.data_out;
    
    models = cell(numel(data_out.stability),1);
SLOPES = nan(numel(data_out.stability),2);
INTERC = SLOPES;
[~,sn,~]=fileparts(files(iF).name);
SESSIONS= repmat({sn},numel(data_out.stability),1);
CLUSTERS = data_out.goodCells';
for iC=1:numel(data_out.stability)
    if startsWith(data_out.region{iC},'RS')
        data_out.region{iC}='RSP';
    end
    if isempty(data_out.region{iC})
        data_out.region{iC}='NULL';
    end
    
    if data_out.stability(iC)>0.4
        x=data_out.allDelays(iC,:,3);
        y=data_out.allDelays(iC,:,1);
        mdl = fitlm(x(:),y(:),'linear');
        INTERC(iC,:)=mdl.Coefficients{1,[1 4]};
        SLOPES(iC,:)=mdl.Coefficients{2,[1 4]};

        
        models{iC}=mdl;
%         plot(x,y,'ro')
%         hold on
%         plot(mdl)
%         pause
%         clf
    else
        models{iC}=nan(1);
    end
end

regs = unique(data_out.region);
for ii=1:numel(regs)
    try
    ridx = strcmp(data_out.region,regs{ii});
    if length(regs{ii})>4
        regs{ii}=regs{ii}(1:4);
    end
    if isfield(ALL_SLOPES,regs{ii})
    ALL_SLOPES.(regs{ii}) = cat(1,ALL_SLOPES.(regs{ii}),SLOPES(ridx,:));
    ALL_SESSIONS.(regs{ii}) = cat(1,ALL_SESSIONS.(regs{ii}),SESSIONS(ridx));
    ALL_CLUSTERS.(regs{ii}) = cat(1,ALL_CLUSTERS.(regs{ii}),CLUSTERS(ridx));
    else
        ALL_SLOPES.(regs{ii})=SLOPES(ridx,:);
        ALL_SESSIONS.(regs{ii}) = SESSIONS(ridx);
        ALL_CLUSTERS.(regs{ii}) = CLUSTERS(ridx);
    end
    catch ME
        disp(ME.message)
    end
end
   
end
%%
regions = {'VISp','MEC'};
figure
hold on
for iR=1:numel(regions)
    dat = ALL_SLOPES.(regions{iR});
    idx = dat(:,2)<0.05;
    histogram(dat(idx,1),'Normalization','probability')
end

legend(regions)
xlabel('Slope')
%%
regions = {'MEC'};
edges = [-0.5:0.05:0.5];
for iR=1:numel(regions)
    dat = ALL_SLOPES.(regions{iR});
    [sess,~,sid]=unique(ALL_SESSIONS.(regions{iR}));
    S_H=zeros(numel(sess),20);
    for iS=1:numel(sess)
    idx = dat(:,2)<0.05 & sid ==iS;
    
    
    [N,edges] = histcounts(dat(idx,1),edges, 'Normalization', 'probability');
    S_H(iS,:)=N;
    end
    %histogram(dat(idx,1),'Normalization','probability')
end
figure
[a,sid]=sort(mean(S_H(:,1:10),2));
imagesc(S_H(sid,:),[0 0.3])
set(gca,'XTick',[1:20],'XTickLabel',.5*edges(1:end-1)+.5*edges(2:end),'XTickLabelRotation',45)


%%
reg ='VISp';
tmp = ALL_SLOPES.(reg)(:,1);

[ma,sid]=sort(abs(tmp),'descend','MissingPlacement','last');
root = 'F:\NP_Data';
ii=17;
while ii<18
    pval = ALL_SLOPES.(reg)(sid(ii),2);
    if pval<0.05
    sn = ALL_SESSIONS.(reg){sid(ii)};
    cid = ALL_CLUSTERS.(reg)(sid(ii));
    data = load(fullfile(root,sn));
    data.posWindow = [300 390];
    
    plotRasterSpeedSort(data,params,[],cid,[5:20])
    xlabel(sprintf('%.2f',tmp(sid(ii))))
    drawnow
    ii=ii+1;
    end
end