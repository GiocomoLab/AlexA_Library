%%
files = dir('F:/NP_DATA/np*dark*');
binsize=10;
DEPTH = [];
MAXF=[];
PI=[];
for iF=1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    if ~ismember('anatomy',fieldnames(data))
        continue
    end
    
    speed = diff(data.posx);

speed(speed<-10) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = [0;speed];
distance = cumsum(speed);
distance(distance<0)=0;
edges=[0:binsize:(max(distance)+binsize)];

distance_trial = ones(size(data.post));




    idxVR=true(size(distance_trial));
    t_time=data.post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=data.sp.st<stop & data.sp.st>=start;
    [spM, dT]=getSpikeMatPosition(data.sp.st(idxNP),data.sp.clu(idxNP),distance,data.post(idxVR),'edges',edges,'max_clust',max(data.sp.clu)+1);

spatialMap=bsxfun(@rdivide,spM(data.sp.cids+1,:,:),dT);
spatialMap(isnan(spatialMap))=0;

smoothSigma = 20/binsize;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

%spatialMap = conv2(spatialMap,gauss_filter','same');


good_cells = data.sp.cgs==2 & strcmp(data.anatomy.cluster_parent,'MEC')';
depth = data.anatomy.tip_distance(good_cells);
[~,sid]=sort(depth,'descend');
zSpatialMap = zscore(spatialMap(good_cells,:),0,[2]);
[pxx,f]=pwelch(zSpatialMap',[],[],[],1/binsize);
[aa,ii]=max(pxx,[],1);

MAXF=cat(1,MAXF,f(ii));
DEPTH = cat(1,DEPTH,depth-data.anatomy.z2);

lags = 1000;
nUnits = size(zSpatialMap,1);
ACG = zeros(nUnits,lags+1);

for ii=1:nUnits
    tmp=xcorr(zSpatialMap(ii,:),1000,'coeff');
    ACG(ii,:)=tmp(1001:end);
end
fig = figure();
[~,sid]=sort(depth,'descend');
imagesc(ACG(sid,1:50),[0 0.5])
hold on
for ii=1:nUnits

    [p,ip]=findpeaks(ACG(sid(ii),1:50),'SortStr','descend');
    if ~isempty(ip)
    PI(end+1)=ip(1);
            plot(ip(1),ii,'kx')

    else
        PI(end+1)=nan;
    end
end
pause
close(fig)

end