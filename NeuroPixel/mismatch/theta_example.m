%% works well with  AA_200123_5_200206_mismatch_1.mat
mm_resp = mean(count_vec(:, 110:130),2)-mean(count_vec(:,80:99),2);
[~,sid]=sort(mm_resp,'descend');

figure
plot(count_vec(sid(1:3),:)')

plot(mean(count_vec(sid(end-10:end),:)))
%%
figure
for ii=1:100
cluID = good_cells(sid(ii));
idx = trial_vec(:,2)==cluID;

tvec = opt.time_bins(1:end-1)*.5+opt.time_bins(2:end)*.5;
        [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
        
        subplot(2,1,1)
        plot(tvec,spike_count)
subplot(2,1,2)
       scatter(trial_vec(idx,1),trial_vec(idx,3),1,'k.')
       pause
       cla
end