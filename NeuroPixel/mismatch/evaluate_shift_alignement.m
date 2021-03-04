%%
mf = dir('/Users/attialex/temp/*.mat');
opt = load_mismatch_opt;
for iF=1:numel(mf)

mm_data = load(fullfile(mf(iF).folder,mf(iF).name));
model_data = load(fullfile('/Users/attialex/model_out/',mf(iF).name));


trial_vec_shifted = mm_data.trial_vec;
random_vec_shifted = mm_data.trial_vec_random;
nT=numel(unique(mm_data.trial_vec(:,3)));

for iT=1:nT
    idx = trial_vec_shifted(:,3)==iT;
    trial_vec_shifted(idx,1) = trial_vec_shifted(idx,1)-model_data.shift_frac(iT);
end

for iT=1:size(model_data.shift_frac_random,2)
    idx = random_vec_shifted(:,3)==iT;
    random_vec_shifted(idx,1) = random_vec_shifted(idx,1)-model_data.shift_frac_random(iT);
end

MM = zeros(numel(mm_data.good_cells),numel(opt.time_bins)-1);
MM_random = MM;
for iC=1:numel(mm_data.good_cells)
    idx = trial_vec_shifted(:,2)==mm_data.good_cells(iC);
    [spike_count]=histcounts(trial_vec_shifted(idx,1),opt.time_bins);
    MM(iC,:)=spike_count;
    idx = random_vec_shifted(:,2)==mm_data.good_cells(iC);
    [spike_count]=histcounts(random_vec_shifted(idx,1),opt.time_bins);
    MM_random(iC,:)=spike_count;
end


FR=mm_data.firing_rate;
MM_ms = MM-mean(MM(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_orig = mm_data.count_vec-mean(mm_data.count_vec(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
MM_r=MM_random-mean(MM_random(:,opt.time_bins>=-.5 & opt.time_bins<0),2);
mm_resp = mean(MM_ms(:,opt.time_bins>0.1 & opt.time_bins<0.6),2);
mm_resp= mm_resp./FR;

valid_idx = FR>1;


figure
plot(opt.time_vecs,mean(MM_ms(valid_idx,:)./FR(valid_idx)))
hold on
%plot(opt.time_vecs,mean(MM_orig(valid_idx,:)./FR(valid_idx)))
plot(opt.time_vecs,mean(MM_r(valid_idx,:)./FR(valid_idx)));
legend({'MM','Random'})

[~,sid]=sort(mm_resp(valid_idx),'descend');
MM_valid = MM_ms(valid_idx,:);
figure
imagesc(MM_valid(sid,:)./FR(valid_idx),[-2 2])
drawnow
end