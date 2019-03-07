correlation_All=zeros(size(spatialMap,3),size(spatialMap,3),size(spatialMap,1));
diagAll=zeros(size(spatialMap,1),size(spatialMap,3)-1);
for iC=1:size(spatialMap,1)
    tmp=corr(squeeze(spatialMap(iC,2:end,:)));
    correlation_All(:,:,iC)=tmp;
    diagAll(iC,:)=diag(tmp,1);
end
%%
blocks=zeros(size(correlation_All));
templates={};
for iC=1:size(correlation_All,3)
    block_i=blockify(squeeze(correlation_All(:,:,iC)),.3,4);
    %block_i=blockify(block_i,.3,2);
    %block_i=blockify(block_i,.3,2);
    cc=bwconncomp(block_i,4);
    for iB=1:cc.NumObjects
        
        [irow,icol]=ind2sub(size(block_i),cc.PixelIdxList{iB});
        mir=min(irow);
        mar=max(irow);
        mic=min(icol);
        mac=max(icol);
        block_i(mic:mac,mir:mar)=1;
        %block_i(irow,icol)=1;    
    end
    blocks(:,:,iC)=double(block_i);
end

%view_stack(cat(2,blocks,correlation_All))

%%
% 
% figure
% for cellIDX=1
% clf
% block=blocks(:,:,cellIDX);
% subplot(1,2,1)
% imagesc(correlation_All(:,:,cellIDX),[0 .5])
% cc=bwconncomp(block,4);
% templates=zeros(40,cc.NumObjects);
% for iB=1:cc.NumObjects
%     [irow,icol]=ind2sub(size(block),cc.PixelIdxList{iB});
%         mir=min(irow);
%         mar=max(irow);
%     templates(:,iB)=squeeze(mean(spatialMap(cellIDX,:,mir:mar),3));
% end
% if cc.NumObjects>1
% y=pdist(templates','correlation');
% Z1=linkage(y,'average');
% subplot(1,2,2)
% dendrogram(Z1)
% end
% saveas(gcf,sprintf('C:/tmp/cell_%d.png',cellIDX))
% end
%%

figure
for cellIDX=1:size(spatialMap,1)
clf
block=blocks(:,:,cellIDX);
cc=bwconncomp(block,4);
templates=zeros(40,cc.NumObjects);
for iB=1:cc.NumObjects
    [irow,icol]=ind2sub(size(block),cc.PixelIdxList{iB});
        mir=min(irow);
        mar=max(irow);
    templates(:,iB)=squeeze(mean(spatialMap(cellIDX,:,mir:mar),3));
end

subplot(2,3,1)
imagesc(squeeze(spatialMap(cellIDX,:,:)))

subplot(2,3,2)
imagesc(correlation_All(:,:,cellIDX),[0 .5])
subplot(2,3,3)
imagesc(blocks(:,:,cellIDX))
xlabel(sprintf('found %d blocks',cc.NumObjects));

if cc.NumObjects>1
y=pdist(templates','correlation');
Z1=linkage(y,'average');
id = cluster(Z1,'cutoff',1);
n_new=length(unique(id));
merged = zeros(size(templates,1),n_new);
for ii=1:n_new
    merged(:,ii)=mean(templates(:,id==ii),2);
end
if n_new>2
sm=squeeze(spatialMap(cellIDX,:,:));

[k,d]=dsearchn(merged',sm');
[a,b]=sort(k);
transitions = find(diff(a)>0);
subplot(2,3,4)
imagesc(sm(:,b))
xdata=[transitions, transitions];
line(xdata',[0 size(sm,1)],'Color','red')
subplot(2,3,5)
imagesc(correlation_All(b,b,cellIDX),[0 .5]);
line(xdata',[0 length(b)],'Color','red')
xlabel(sprintf('ended up with %d unique blocks',n_new))
else
    if cc.NumObjects>1

    y=pdist(templates','correlation');
Z1=linkage(y,'average');
subplot(2,3,4)
dendrogram(Z1)
xlabel('only 1 block after merging')
    end
end

end
saveas(gcf,sprintf('C:/tmp/cell_%d.png',cellIDX))

end