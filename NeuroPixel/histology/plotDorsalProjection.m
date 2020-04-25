bregma=allenCCFbregma();
figure
hold on
for iA=1:numel(AllProbePoints)
    currentProbes = AllProbePoints{iA};
    for iP = 1:numel(currentProbes)
        m = currentProbes{iP};
        m=m-bregma;
        plot(m(3),m(1),'o')
    end
end
set(gca,'YDir','reverse')

%%
pp = av(:,1:20,:);
pp=median(pp,2);

figure
imagesc(squeeze(pp))
%%
xdim = size(av,1);
ydim = size(av,3);
Updated = false(xdim,ydim);
projection = zeros(xdim,ydim);


for iZ = 1:size(av,2)
    plane = squeeze(av(:,iZ,:));
    notOne = plane~=1;
    toUpdate = notOne & ~Updated;
    projection(toUpdate)=plane(toUpdate);
    Updated = Updated | notOne;
end
%%
figure
imagesc(projection)
hold on


for iA=1:numel(AllProbePoints)
    currentProbes = AllProbePoints{iA};
    for iP = 1:numel(currentProbes)
        m = currentProbes{iP};
        %m=m-bregma;
        plot(m(3),m(1),'rx')
    end
end

%%
figure('Position',[680   306   561   672],'Renderer','Painters')
outlines = zeros(size(projection));
uA= unique(projection);
for iA=2:numel(uA)
BW = projection == uA(iA);
[B,L] = bwboundaries(BW,'noholes');
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1),'Color', [.2 .2 .2], 'LineWidth', 1)
end
end


for iA=1:numel(AllProbePoints)
    currentProbes = AllProbePoints{iA};
    for iP = 1:numel(currentProbes)
        m = currentProbes{iP};
        %m=m-bregma;
        plot(m(3),m(1),'ro')
    end
end

set(gca,'YDir','reverse')

for iP = 1:numel(AllPointsMEC)
    plot(AllPointsMEC{iP}(2),AllPointsMEC{iP}(1),'ro')
end
    