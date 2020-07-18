%% plot pca of pb resp and mm resp, pca based on concatenated matrix
opt = load_mismatch_opt;
opt.time_bins =-2:0.02:3;
opt.time_vecs = opt.time_bins(1:end-1)*0.5+opt.time_bins(2:end)*0.5;

tidx = opt.time_vecs>-.4 & opt.time_vecs<1;
pre_window = opt.time_vecs<-0.1 & opt.time_vecs>-.5;
n_samples_per_cond = nnz(tidx);
uS = unique(SID);
mm_resp = MM(:,tidx)-mean(MM(:,pre_window),2);
pb_resp = PB(:,tidx)-mean(PB(:,pre_window),2);

smooth_method = 'gaussian';
smooth_window = 9;
mm_resp = smoothdata(mm_resp,2,smooth_method,smooth_window);
pb_resp = smoothdata(pb_resp,2,smooth_method,smooth_window);

mm_resp = zscore(mm_resp,[],2);
pb_resp = zscore(pb_resp,[],2);
for iS=1:numel(uS)
    sidx = SID==uS(iS);
    pb_resp_this = pb_resp(sidx,:);
r_idx = randperm(size(pb_resp_this,1));
%pb_resp_this=pb_resp_this(r_idx,:);


[coeff,score,latent,tsquared,explained,mu] = pca(mm_resp(sidx,:)');

tmp=mm_resp(sidx,:)';
tmp=tmp-mu;
sc1=tmp*coeff(:,1:3);
p = ALLMEC(iS,idx_u)-pca_data.pca_mu;
    c = p*pca_data.pca_coeff(:,1:3);
 %color = cat(1,hsv(n_samples_per_cond));
   
% figure('Color','White')
% scatter3(score(:,1),score(:,2),score(:,3),35,color,'.')
% hold on
% patch([score(:,1); nan],[score(:,2); nan],[score(:,3); nan],[opt.time_vecs(tidx) nan],'FaceColor','none','EdgeColor','interp')
% colorbar
% colormap(hsv)
% view(3)
% 
data_matrix = cat(1,mm_resp(sidx,:)',pb_resp_this');
color = cat(1,hsv(n_samples_per_cond),hsv(n_samples_per_cond));

[coeff,score,latent,tsquared,explained,mu] = pca(data_matrix);
score(n_samples_per_cond,:)=nan;
figure('Color','White')
scatter3(score(:,1),score(:,2),score(:,3),45,color,'.')
hold on
patch([score(:,1); nan],[score(:,2); nan],[score(:,3); nan],[opt.time_vecs(tidx) opt.time_vecs(tidx) nan],'FaceColor','none','EdgeColor','interp','LineWidth',2)
colormap(hsv)
end
%% pb resp projected onto MM pca
tidx = opt.time_vecs>-.4 & opt.time_vecs<1;
pre_window = opt.time_vecs<-0.1 & opt.time_vecs>-.5;
n_samples_per_cond = nnz(tidx);
uS = unique(SID);
mm_resp = MM(:,tidx)-mean(MM(:,pre_window),2);
pb_resp = PB(:,tidx)-mean(PB(:,pre_window),2);

smooth_method = 'gaussian';
smooth_window = 9;
mm_resp = smoothdata(mm_resp,2,smooth_method,smooth_window);
pb_resp = smoothdata(pb_resp,2,smooth_method,smooth_window);

mm_resp = zscore(mm_resp,[],2);
pb_resp = zscore(pb_resp,[],2);

for iS=1:numel(uS)
    figure('Color','White')

    sidx = SID==uS(iS);
    pb_resp_this = pb_resp(sidx,:);
r_idx = randperm(size(pb_resp_this,1));
%pb_resp_this=pb_resp_this(r_idx,:);


[coeff,score,latent,tsquared,explained,mu] = pca(mm_resp(sidx,:)');

tmp= pb_resp_this'-mu;
sc1=tmp*coeff(:,1:3);
tmp= pb_resp_this(r_idx,:)'-mu;
sc2=tmp*coeff(:,1:3);
 color_mm = cat(1,hsv(n_samples_per_cond+1));
 color_mm = cbrewer('seq','Reds',n_samples_per_cond+1);
 color_mm = reshape(color_mm,[1,size(color_mm,1),size(color_mm,2)]);
 
 color_pb = cbrewer('seq','Greens',n_samples_per_cond+1);
 color_pb = reshape(color_pb,[1,n_samples_per_cond+1,3]);
  color_pbr = cbrewer('seq','Blues',n_samples_per_cond+1);
 color_pbr = reshape(color_pbr,[1,n_samples_per_cond+1,3]);


%scatter3(score(:,1),score(:,2),score(:,3),35,color,'.')
patch([score(:,1);nan],[score(:,2); nan],[score(:,3); nan],color_mm,'FaceColor','none','EdgeColor','interp','LineWidth',2)
set(gca, 'Projection','perspective')
view(3)
hold on
patch([sc1(:,1); nan],[sc1(:,2); nan],[sc1(:,3); nan],color_pb,'FaceColor','none','EdgeColor','interp','LineWidth',2)
patch([sc2(:,1); nan],[sc2(:,2); nan],[sc2(:,3); nan],color_pbr,'FaceColor','none','EdgeColor','interp','LineWidth',2)

%pause
%clf
end
%%
%% pb resp projected onto MM pca, across all
    figure('Color','White')

tidx = opt.time_vecs>-.4 & opt.time_vecs<1;
pre_window = opt.time_vecs<-0.1 & opt.time_vecs>-.5;
n_samples_per_cond = nnz(tidx);
uS = unique(SID);
mm_resp = MM(:,tidx)-mean(MM(:,pre_window),2);
pb_resp = PB(:,tidx)-mean(PB(:,pre_window),2);

smooth_method = 'gaussian';
smooth_window = 9;
mm_resp = smoothdata(mm_resp,2,smooth_method,smooth_window);
pb_resp = smoothdata(pb_resp,2,smooth_method,smooth_window);

mm_resp = zscore(mm_resp,[],2);
pb_resp = zscore(pb_resp,[],2);


r_idx = randperm(size(pb_resp,1));



[coeff,score,latent,tsquared,explained,mu] = pca(mm_resp');

tmp= pb_resp'-mu;
sc1=tmp*coeff(:,1:3);
tmp= pb_resp(r_idx,:)'-mu;
sc2=tmp*coeff(:,1:3);
 color_mm = cat(1,hsv(n_samples_per_cond+1));
 color_mm = cbrewer('seq','Reds',n_samples_per_cond+1);
 color_mm = reshape(color_mm,[1,size(color_mm,1),size(color_mm,2)]);
 
 color_pb = cbrewer('seq','Greens',n_samples_per_cond+1);
 color_pb = reshape(color_pb,[1,n_samples_per_cond+1,3]);
  color_pbr = cbrewer('seq','Blues',n_samples_per_cond+1);
 color_pbr = reshape(color_pbr,[1,n_samples_per_cond+1,3]);


%scatter3(score(:,1),score(:,2),score(:,3),35,color,'.')
patch([score(:,1);nan],[score(:,2); nan],[score(:,3); nan],color_mm,'FaceColor','none','EdgeColor','interp','LineWidth',2)
set(gca, 'Projection','perspective')
view(3)
hold on
patch([sc1(:,1); nan],[sc1(:,2); nan],[sc1(:,3); nan],color_pb,'FaceColor','none','EdgeColor','interp','LineWidth',2)
patch([sc2(:,1); nan],[sc2(:,2); nan],[sc2(:,3); nan],color_pbr,'FaceColor','none','EdgeColor','interp','LineWidth',2)

%pause
%clf
