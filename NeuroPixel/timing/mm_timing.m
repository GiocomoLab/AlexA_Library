avgMM =aggregateData.avgMM;
mmresp = mean(avgMM(:,220:250),2)-mean(avgMM(:,155:185),2);
tvec = (-200:200)/50;
[~,sidx]=sort(mmresp);
figure
subplot(2,1,1)
plot(tvec,avgMM(sidx(1:3),:)')
grid on
xlim([-0.5 1])
subplot(2,1,2)
[~,sidx]=sort(mmresp,'descend');

plot(tvec,avgMM(sidx(1:3),:)')
grid on
xlim([-0.5 1])
%%
    iSite = 1;
    speed_t=0.05;
    nTot = 0;
    n_sess= size(aggregateData.MM_snps{iSite},1);
    thisSession_idx = nTot+1:nTot+n_sess;
    avgMM = aggregateData.avgMM(thisSession_idx,:);
    mmresp = mean(avgMM(:,220:250),2)-mean(avgMM(:,155:185),2);
    
    [~,sidx]=sort(mmresp,'descend');

    win=[-4 4];
    spike_mat = aggregateData.MM_snps{iSite};
    binned_array=squeeze(spike_mat(sidx(1),:,:));
    bins=win(1):0.001:win(2);
    
    % set the new rasters
%     [tr,b] = find(binned_array);
%     [rasterX,yy] = rasterize(bins(b));
%     rasterY = yy*1+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
%     figure;
%     plot(rasterX,rasterY)
    
 % 1 filter for true MM events (based on speed)
 mm_run = aggregateBeh.MMAllRun{iSite};
run_bi = mm_run>speed_t;
valid_idx = all(run_bi(:,175:225),2);

binned_array = spike_mat(:,valid_idx,:);
%turn into firing rate for 20ms chunks
stepsize = 20;
rate_mat = zeros(size(binned_array,1),nnz(valid_idx),floor(8001/20));
cntr = 0;
for iS = 1:stepsize:(8001-stepsize)
    cntr = cntr+1;
    idx = iS:iS+stepsize;
    n_spikes =  sum(binned_array(:,:,idx),3);
    rate_mat(:,:,cntr)=n_spikes;
end
sig_val = zeros(size(rate_mat,1),size(rate_mat,3));

for iC = 1:size(rate_mat,1)
for iS=1:size(sig_val,2)

pre_dist = reshape(rate_mat(iC,:,(2000/stepsize:4000/stepsize)),1,[]);

    [h,p]=ttest2(pre_dist,squeeze(rate_mat(iC,:,iS)));
    sig_val(iC,iS)=p;
end
end
%%
win=4000/stepsize:5000/stepsize;
%win = 1:size(rate_mat,3);
crossings = nan(size(sig_val,1),1);
for iC= 1:size(rate_mat,1)
    a=strfind(sig_val(iC,win)<0.05,[0 1 1])+1;
    if ~isempty(a)
        crossings(iC)=a(1);
    end
end
figure
tvec = linspace(-4,4,size(rate_mat,3));
tmp = mean(squeeze(mean(rate_mat(sidx(1:20),:,:),2)));
plot(tvec,tmp)
hold on
ff=nanmedian(crossings(sidx(1:20)));
ff=round(ff)+win(1)-1;
plot(tvec(ff),tmp(ff),'r*')
%lot(
% figure
% imagesc(sig_val<0.05)
% hold on
% plot(crossings+win(1)-1,1:numel(crossings),'r*')

 % for each time bin take pre distribution as null distribution, then run
 % ttest2 for each time bin.
 
 %suggestion by kiah: shuffle data and get a null distribution like that.