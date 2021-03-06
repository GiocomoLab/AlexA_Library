

%% load data

%%
try
    if isfield(data.anatomy,'parent_shifted')
        region = data.anatomy.parent_shifted;
    else
        region = data.anatomy.cluster_parent;
    end
catch
    disp('no anatomy')
    return
end

%%
u_regs = unique(region);
c = arrayfun(@(x)sum(strcmp(region,x)), unique(region), 'Uniform', false);
counts = cell2mat(c);

errors={};
location_real={};
location_decoded={};
n_cells{iR}={};
gains ={};
%% select cells
%good_cells = data.sp.cids(data.sp.cgs==2 & ismember(region,'VISp'));
%good_cells = data.sp.cids(ismember(data.anatomy.parent_shifted,'VISp'));
for iR=1:numel(u_regs)
    current_reg = u_regs{iR};
good_cells = data.sp.cids(data.sp.cgs>=1 & ismember(region,current_reg));
if numel(good_cells)<50
    continue
end
% calc spatial Firing rate map for these cells for the selected trial range

idxClu = ismember(data.sp.clu,(good_cells));

idxVR=data.trial<=20;
t_time=data.post(idxVR);
start=min(t_time);
stop=max(t_time);
idxNP=data.sp.st<stop & data.sp.st>=start;

[sp,occ,good_cells,track_edges]=getSpikeMatPosition2(data.sp.st(idxClu&idxNP),data.sp.clu(idxClu&idxNP),data.posx(idxVR),data.post(idxVR));
sp=reshape(sp,[size(sp,1) 1 size(sp,2)]);

% calculate firing rate by time

tBin = 0.25;
nBins = ceil(max(data.sp.st)/tBin);
edges=0:tBin:nBins*tBin;
tvec=edges(1:end-1)+edges(2:end);
tvec=tvec/2;
frMat = zeros(nBins+1,numel(good_cells));
filt = gausswin(5);
filt = filt/sum(filt);
for iC=1:numel(good_cells)
    frMat(:,iC)=conv(histc(data.sp.st(data.sp.clu==good_cells(iC)),edges),filt,'same');
end

% calculate posterior
withHistory = false;
post = decode_calcBayesPost(frMat', sp, occ',tBin,withHistory);
% decode position for selected trials
test_trials = [21:24];
idxVR=ismember(data.trial,test_trials);
    t_time=data.post(idxVR);
    start=min(t_time);
    stop=max(t_time);
idxPost = tvec<stop & tvec>start;
figure
subplot(1,2,1)
imagesc(flipud(log(squeeze(post(:,:,idxPost)))))
subplot(1,2,2)
tmp = data.posx(idxVR);
tmp(tmp<0)=0;
tmp(tmp>=400)=400;
aa=discretize(tmp,track_edges);

%plot(data.post(idxVR),data.posx(idxVR))
plot(data.post(idxVR),aa)
hold on
[~,maxi]=max(post(:,:,idxPost),[],1);
plot(tvec(idxPost),(squeeze(maxi)))
track_centers = track_edges(1:end-1)+track_edges(2:end);
track_centers = track_centers /2;
loc_downsampled = interp1(data.post(idxVR),track_centers(aa),tvec(idxPost));
error = (loc_downsampled-track_centers(maxi));

errors{iR}=error;
location_real{iR}=loc_downsampled;
location_decoded{iR}=track_centers(maxi);
n_cells{iR}=numel(good_cells);
gains{iR}=data.trial_gain(test_trials);
end