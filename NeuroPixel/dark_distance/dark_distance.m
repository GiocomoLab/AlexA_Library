spike_id=sp.clu==140;
spike_t = sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,post);
figure
scatter(posx(spike_idx),trial(spike_idx),2)

%%
dp2=diff(posx);
teleport = find(dp2<-1);
dp2(teleport)=0.5*dp2(teleport-1)+0.5*dp2(teleport+1);
distance = [0; cumsum(dp2)];

%%
%spike_id=ismember(sp.clu,sp.cids(sp.cgs==2));

spike_t = sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,post);
spike_distance = distance(spike_idx);


[CGR,b]=CCG(spike_distance,double(sp.clu(spike_id)),'binSize',[2],'duration',[10000]);
figure
plot(b,squeeze(CGR(:,140,140)))
binedges = 0:2:max(distance)+2;
time_per_bin = histcounts(distance, binedges);
time_per_bin = time_per_bin * params.TimeBin;
[aa,bb]=histcounts(spike_distance,binedges);
firing_rate = aa./time_per_bin;
if sum(isnan(firing_rate))>0
    firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
end
smoothSigma = params.SmoothSigmaFR/params.SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
fr_smooth = conv(firing_rate,gauss_filter,'same');
[cc,dd]=xcorr(fr_smooth,10000);
figure
plot(b,squeeze(CGR(:,140,140))/max(squeeze(CGR(:,140,140))))
hold on
plot(dd*2,cc/max(cc),'r')
%%
[aa,bb]=histcounts(spike_distance,0:2:max(distance)+4);

[pxx,f]=pwelch(aa,[],[],[],1/2);
figure
plot(f,pxx)
%%
good_cells = sp.cids(sp.cgs==2);
ACG = zeros(201,numel(good_cells));
for k=1:length(good_cells)
    cluid = good_cells(k);
    spike_id = sp.clu==cluid;
    spike_t = sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,post);
spike_distance = distance(spike_idx);


[CGR,b]=CCG(spike_distance,double(sp.clu(spike_id)),'binSize',[10],'duration',[2000]);
    if nnz(spike_id)>10
    ACG(:,k)=squeeze(CGR(:,end,end));
    end
end
%%
[CGR,b]=CCG(sp.st(spike_id),double(sp.clu(spike_id)),'binSize',[0.1],'duration',[100]);
figure
plot(b,squeeze(CGR(:,140,140)))
%%
% shifts = randi([-1900 1900],50,1);
% ACG = zeros(501,numel(shifts));
% spike_id = sp.clu==140;
%     spike_t = sp.st(spike_id);
% [~,~,spike_idx] = histcounts(spike_t,post);
% spike_distance = distance(spike_idx);
% for k=1:length(shifts)
%     tmp_dist = circshift(spike_distance,shifts(k),1);
%     [CGR,b]=CCG(tmp_dist,double(sp.clu(spike_id)),'binSize',[20],'duration',[10000]);
%     ACG(:,k)=squeeze(CGR(:,end,end));
% end

%%
good_cells = sp.cids(sp.cgs==2);
PXX = zeros(4097,numel(good_cells));
time_per_bin = histcounts(distance, binedges);
binedges = 0:2:max(distance)+2;
ACG=zeros(numel(good_cells),10001);
time_per_bin = time_per_bin * params.TimeBin;
for k=1:length(good_cells)
    cluid = good_cells(k);
    spike_id = sp.clu==cluid;
    spike_t = sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,post);
spike_distance = distance(spike_idx);


[aa,bb]=histcounts(spike_distance,binedges);
firing_rate = aa./time_per_bin;
if sum(isnan(firing_rate))>0
    firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
end
smoothSigma = 3;%params.SmoothSigmaFR/params.SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
fr_smooth = conv(firing_rate,gauss_filter,'same');



[pxx,f]=pwelch(fr_smooth,[],[],[],1/2);
PXX(:,k)=pxx;
[cc,dd]=xcorr(fr_smooth,5000);
ACG(k,:)=cc;
end
%%
idx = f>8.5449e-04;
idx(end)=false;
start = find(idx,1);
start=start-1;

[A,B]=sort(PXX(idx,:),'descend');
PXXtmp = PXX(idx,:);
%for ic=1:384; PXXtmp(:,ic)=PXXtmp(B(:,ic),ic); end;
diff_left = zeros(size(A));
diff_right = zeros(size(A));
for ic=1:384;
diff_left = A(:,ic)-PXX(B(:,ic)+start-1,ic);
diff_right = A(:,ic)-PXX(B(:,ic)+start+1,ic);
end

pp = 1./f;
figure
histogram(pp(ii),0:10:1000)

%%
figure

for k=1:length(good_cells)
    cluid = good_cells(k); 
    spike_id = sp.clu==cluid;
    spike_t = sp.st(spike_id);
    [~,~,spike_idx] = histcounts(spike_t,post);
    spike_distance = distance(spike_idx);
    phase = pp(ii(k));
    norm_dist = mod(spike_distance,phase);
    norm_trial = floor(spike_distance/phase);
    scatter(norm_dist,norm_trial)
    axis tight
    pause
    clf
end
%%
figure
idx = f>8.5449e-04;
idx(end)=false;
start = find(idx,1);
start=start-1;
locs = zeros(2,size(PXX,2));
for k=1:length(good_cells)
    tmp = PXX(:,k);
    med = median(tmp);
    [a,b]=sort(PXX(idx,k),'descend');
    diff_left = a-PXX(b+start-1,k);
    diff_right = a-PXX(b+start+1,k);
    plot(f,PXX(:,k))
    hold on
    pot_idx = find(diff_left>0 & diff_right > 0,2);
    maxloc = b(pot_idx(1))+start;
    sec_loc = b(pot_idx(2))+start;
    plot(f(maxloc),PXX(maxloc,k),'ro')
    plot(f(sec_loc),PXX(sec_loc,k),'go')
    xlim([-1000 1000])
    axis tight
    pause
    clf
end
%%

idx = f>8.5449e-04;
idx(end)=false;
start = find(idx,1);
start=start-1;
locs = zeros(2,size(PXX,2));
for k=1:length(good_cells)
    tmp = PXX(:,k);
    med = median(tmp);
    [a,b]=sort(PXX(idx,k),'descend');
    diff_left = a-PXX(b+start-1,k);
    diff_right = a-PXX(b+start+1,k);
   
    pot_idx = find(diff_left>0 & diff_right > 0,2);
    maxloc = b(pot_idx(1))+start;
    sec_loc = b(pot_idx(2))+start;
    locs(:,k)=[maxloc,sec_loc];
end

pp = 1./f;
figure
histogram(pp(locs(1,:)),0:10:1000)
%%

figure
idx = (dd*2)>25 & (dd*2)<1000;
start = find(idx,1);
start=start-1;
for k=1:length(good_cells)
    tmp = ACG(k,:);
    med = median(tmp);
    [a,b]=sort(ACG(k,idx),'descend');
    diff_left = a-ACG(k,b+start-1);
    diff_right = a-ACG(k,b+start+1);
    plot(dd*2,ACG(k,:))
    hold on
    pot_idx = find(diff_left>0 & diff_right > 0,2);
    maxloc = b(pot_idx(1))+start;
    sec_loc = b(pot_idx(2))+start;
    plot(dd(maxloc)*2,ACG(k,maxloc),'ro')
    plot(dd(sec_loc)*2,ACG(k,sec_loc),'go')
    xlim([-1000 1000])
    axis tight
    pause
    clf
end
