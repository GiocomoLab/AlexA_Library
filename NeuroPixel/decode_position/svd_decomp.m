%data = load('F:\NP_DATA\AA2_190810_contrast_1.mat');
data = load('F:\NP_DATA\npF4_1025_gaincontrast_2.mat')
%data=load('F:\Dropbox\pretty_rasters_whole_session_combined\npF4_1025_gaincontrast_2_all_rasters_combined.png');
%%
load('Z:\giocomo\attialex\cellInfoGain.mat')
%%
session_name = 'npF4_1025_gaincontrast_2';
local_stab = cell_info.LocalStability;
numblocks = floor(32/4);
local_stab = mean(local_stab(:,1:numblocks),2);
this_sess = strcmp(cell_info.Session,session_name);
this_cells = cell_info.CellID(this_sess);
this_stab = local_stab(this_sess);
this_reg = cell_info.BrainRegion(this_sess);
this_MEC = strcmp(this_reg,'MEC');
stable_cells=this_cells(this_stab>=0.3);
%%

trials=[1:30];
%trials = trials(trial_gain == 1 & trial_contrast == 100);
spatialMap=[];
dwell_time=[];
edges=[0:5:400];
edges(1)=-.01;
data.posx(data.posx<0)=0;
data.posx(data.posx>=400)=399.00;
for iT=1:length(trials)
    idxVR=data.trial==trials(iT);
    t_time=data.post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=data.sp.st<stop & data.sp.st>=start;
    [spM, dT]=getSpikeMatPosition(data.sp.st(idxNP),data.sp.clu(idxNP),data.posx(idxVR),data.post(idxVR),'edges',edges,'max_clust',max(data.sp.clu)+1);
    spatialMap=cat(3,spatialMap,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
reg = strcmp(data.anatomy.cluster_parent,'MEC')';
spatialMap=spatialMap(data.sp.cids+1,:,:);
spatialMap=spatialMap(data.sp.cgs==2,:,:);
spatialMap = spatialMap(this_stab>0.0 & this_MEC,:,:);
%normalize by dwell time in each bin
dt=dwell_time';
dt=reshape(dt,[1 size(dt,1),size(dt,2)]);
for ii=1:size(spatialMap,1)
    spatialMap(ii,:,:)=spatialMap(ii,:,:)./dt;
end
spatialMap(isnan(spatialMap))=0;
%%

X=zeros(size(spatialMap,1),size(spatialMap,2)*size(spatialMap,3));
filt = gausswin(15);
filt = filt/sum(filt);
for iC = 1:size(spatialMap,1)
    tmp = spatialMap(iC,:,:);
    
    X(iC,:)=tmp(:);
end
X(isnan(X))=0;
X=conv2(X,filt','same');
X=X-mean(X,2);
X=normc(X);
%X=X-mean(X,2);
%X=bsxfun(@rdivide,X,std(X,[],2));
%X=bsxfun(@minus,X,mean(X,2));

%%

[U,S,V]=svd(X);

%%
Xp=U(:,1:3)'*X;
filt=gausswin(7);
filt = repmat(filt',[3,1]);
Xp=conv2(Xp,filt,'same');
%%
cm = winter(4);
cm = cat(1,cm,summer(4));
%cm = zeros(20,3);
%cm(:,3)=1;
figure
hold on
nBins = size(spatialMap,2);
for iT=1:24
    idx = (iT-1)*nBins+1:iT*nBins;
    %plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'.','Color',cm(iT,:))
    plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'-','Color',cm(iT,:))
end


% cm = jet(20);
% cm = zeros(20,3);
% cm(:,1)=1;
% for iT=21:24
%     idx = (iT-1)*nBins+1:iT*nBins;
%     plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'.','MarkerEdgeColor',cm(iT-20,:))
% end
%%
figure
bl_traj = zeros(3,size(spatialMap,2));
for iT=4:20
    idx = (iT-1)*nBins+1:iT*nBins;
    bl_traj=bl_traj+1/20*Xp(:,idx);
end
plot3(bl_traj(1,:),bl_traj(2,:),bl_traj(3,:))
hold on
gain_traj = zeros(3,size(spatialMap,2));
for iT=21:24
    idx = (iT-1)*nBins+1:iT*nBins;
    plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'.','MarkerEdgeColor',cm(iT-20,:))
end

%%
figure
subplot(1,2,1)
imagesc(X)
[~,sid]=sort(U(:,1));
subplot(1,2,2)
imagesc(X(sid,:))
%%
[coeff,score,latent,tsquared] = pca(X');
Xp = coeff(:,1:3)'*X;
filt=gausswin(5);
filt = repmat(filt',[3,1]);
Xp=conv2(Xp,filt,'same');
%%
cm = winter(20);
%cm = zeros(20,3);
%cm(:,3)=1;
figure
hold on
nBins = size(spatialMap,2);
for iT=5:30
    idx = (iT-1)*nBins+1:iT*nBins;
    %plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'.','Color',cm(iT,:))
    plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'-','Color',cm(iT,:))
end


cm = jet(20);
cm = zeros(20,3);
cm(:,1)=1;
for iT=21:24
    idx = (iT-1)*nBins+1:iT*nBins;
    %plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'.','MarkerEdgeColor',cm(iT-20,:))
    plot3(Xp(1,idx),Xp(2,idx),Xp(3,idx),'r-')
end
