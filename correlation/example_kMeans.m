filenames = {'G4/1204_mismatch_1/1204_mismatch_1.mat',...
    'G2/1211_mismatch_1/1211_mismatch_1.mat',...
    'G2/1212_mismatch_1/1212_mismatch_1.mat',...
    'G5/1207_mismatch_1/1207_mismatch_1.mat',...
    'G5/1210_mismatch_1/1210_mismatch_1.mat'
    };


root_dir='F:\';

CGR_all=[];
for iF=4
    load([root_dir filenames{iF}]);
    good_cells = sp.cids(sp.cgs==2);
    idx=ismember(sp.clu,good_cells);
    [CGR,b]=CCG(sp.st(idx),double(sp.clu(idx)),'binSize',[0.001],'duration',[0.5]);
    CGRGood=CGR(:,good_cells,good_cells);
    
    %%
    CGR_cat = reshape(CGRGood,size(CGR,1),[]);
    
    stride = size(CGRGood,2);
    keep=true(stride);
    for ii = 1:stride
        IDX = ii+(ii-1)*stride;
        keep(IDX)=0;
    end
    linear=keep(:);
    CGR_cat=CGR_cat(:,linear);
    CGR_all=cat(1,CGR_all,CGR_cat');
end

%first scale by max
CGR_n=CGR_all-mean(CGR_all,2);%./repmat(max(CGR_all,[],2),1,51);
%then whiten each
CGR_ms=CGR_n-mean(CGR_n);
%%
dist=[];
for kk=1:10
    [IDX,~,sumd]=kmeans(CGR_ms,kk);
    dist(end+1)=sum(sumd);
end

figure;plot(dist)
%%
[IDX,~,sumd]=kmeans(CGR_ms,11);
%%
figure('Name','KMEANS')
hold on
for ii=1:9
    subplot(5,2,ii)
    plot([-.025:0.001:.025],mean(CGR_ms(IDX==ii,:)))
    xlabel('Time [s]')
    grid on
    axis tight
end
%%
[coeff,score,latent,~,vex] = pca(CGR_ms);
%%
figure('Name','PCA');
for ii=1:10
    subplot(5,2,ii)
    plot([-.025:0.001:.025],coeff(:,ii))
    xlabel(sprintf('Time [s], Explained: %.2f',vex(ii)))
    grid on
    axis tight
end
%%
% Y=tsne(CGR_ms);
% %%
% gscatter(Y(:,1),Y(:,2),IDX,eye(3))

%%

Y=tsne(score(:,1:20));
%%
nComp=20;

[IDX,~,sumd]=kmeans(score(:,1:nComp) * coeff(:,1:nComp)',11);

figure
gscatter(Y(:,1),Y(:,2),IDX,eye(3))
%%
reconst = score(:,1:nComp) * coeff(:,1:nComp)';

