data=load('/Volumes/T7/AA_201025_3_TowerTraining_201125_18-27-26.mat');
good_cells = data.sp.cids(data.sp.cgs==2);
opt = load_default_opt;
tmp = diff(data.posx)<-50;
data.trial = [0;cumsum(tmp)];
[corrMat,frMat,shiftMat] = trialCorrMat(good_cells,0:max(data.trial)-1,data,opt);
%%
figure
for ii=1:201
    imagesc(squeeze(frMat(ii,:,1:100)))
    pause
    clf
end

%%
trigs_1=strfind(data.vr_data_resampled.Object==1,[0 0 1 1])+2;
trigs_2=strfind(data.vr_data_resampled.Object==0,[0 0 1 1])+2;
opt = load_mismatch_opt;
[spikeTimes1,~,~,~,count_vec1]=extract_triggered_spikeTimes(data.sp,data.post(trigs_1),'cluIDs',good_cells,'win',opt.extract_win);

[spikeTimes2,~,~,~,count_vec2]=extract_triggered_spikeTimes(data.sp,data.post(trigs_2),'cluIDs',good_cells,'win',opt.extract_win);
%%
figure
for ii=1:201
    plot(opt.time_vecs,count_vec1(ii,:))
    hold on
    plot(opt.time_vecs,count_vec2(ii,:))
    pause
    cla
end
%%
cv1_even = spikeStruct2countVec(cat(1,spikeTimes1{1:2:end}),good_cells,opt);
cv1_odd = spikeStruct2countVec(cat(1,spikeTimes1{2:2:end}),good_cells,opt);
pre_window = opt.time_vecs<-0.1 & opt.time_vecs>-.5;
post_window = opt.time_vecs>0.1 & opt.time_vecs<3;
cv1_resp = mean(cv1_even(:,pre_window),2)-mean(cv1_even(:,post_window),2);
[~,sid]=sort(cv1_resp);
figure
subplot(1,2,1)
resp_1=cv1_odd-mean(cv1_odd(:,pre_window),2);
resp_1 = resp_1./nanmedian(count_vec1,2);
imagesc(smoothdata(resp_1(sid,:),2,'gaussian',13),[-1 1])
subplot(1,2,2)
resp_2 = count_vec2-mean(count_vec2(:,pre_window),2);
resp_2 = resp_2./nanmedian(count_vec2,2);
imagesc(smoothdata(resp_2(sid,:),2,'gaussian',13),[-1 1])
%%
figure
for ii=1:201
    plot(opt.time_vecs,cv1_odd(sid(ii),:),'b')
    hold on
    plot(opt.time_vecs,cv1_even(sid(ii),:),'b--')
    hold on
    plot(opt.time_vecs,count_vec2(sid(ii),:))
    pause
    cla
end

%% pca analysis

opt.time_bins =-2:0.02:3;
opt.time_vecs = opt.time_bins(1:end-1)*0.5+opt.time_bins(2:end)*0.5;

tidx = opt.time_vecs>-1 & opt.time_vecs<1;
pre_window = opt.time_vecs<-0.1 & opt.time_vecs>-.5;
n_samples_per_cond = nnz(tidx);
MM=count_vec1;
PB=count_vec2;
mm_resp = MM(:,tidx);%-mean(MM(:,pre_window),2);
pb_resp = PB(:,tidx);%-mean(PB(:,pre_window),2);

smooth_method = 'gaussian';
smooth_window = 9;
mm_resp = smoothdata(mm_resp,2,smooth_method,smooth_window);
pb_resp = smoothdata(pb_resp,2,smooth_method,smooth_window);

%mm_resp = zscore(mm_resp,[],2);
%pb_resp = zscore(pb_resp,[],2);

data_matrix = cat(1,mm_resp',pb_resp');
color = cat(1,hsv(n_samples_per_cond),hsv(n_samples_per_cond));
color = cat(1,brewermap(n_samples_per_cond,'*RdBu'),brewermap(n_samples_per_cond,'*RdBu'));
[coeff,score,latent,tsquared,explained,mu] = pca(data_matrix);
score(n_samples_per_cond,:)=nan;
figure('Color','White')
scatter3(score(:,1),score(:,2),score(:,3),45,color,'.')
hold on
patch([score(:,1); nan],[score(:,2); nan],[score(:,3); nan],[opt.time_vecs(tidx) opt.time_vecs(tidx) nan],'FaceColor','none','EdgeColor','interp','LineWidth',2)
colormap(brewermap(n_samples_per_cond,'*RdBu'))
