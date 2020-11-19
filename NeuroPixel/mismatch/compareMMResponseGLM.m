
good_cells =data.sp.cids(data.sp.cgs==2);
%glmData = fitGLM_mismatch(data,good_cells,10);
glmData = fitGLM_mismatchNew(data,good_cells,10);

%%
[mm,runO]=preview_MM(data);
frMat = calcTrialFRMat(good_cells,1:max(data.trial)-1,data,load_default_opt); % single trial fr mat

%%
opt = load_mismatch_opt;
mm_resp = nanmean(mm(:,opt.time_vecs> 0.1 & opt.time_vecs<1.0),2)-nanmean(mm(:,opt.time_vecs>-1 & opt.time_vecs<-0.1),2);
[~,sidx] = sort(mm_resp,'descend');
figure
for ii=1:numel(good_cells)
    iC=sidx(ii);
    subplot(1,2,1)
    plot(squeeze(nanmean(frMat(iC,:,:),2)))
    hold on
    if ismember(1,glmData(iC).bestModels)
        plot(glmData(iC).tuning_curves{1})
    end
    title(glmData(iC).bestModels)
    xlabel(glmData(iC).final_pval)
    subplot(1,2,2)
    plot(opt.time_vecs,mm(iC,:))
    yyaxis right
    plot(opt.time_vecs,runO(iC,:))
    pause
    clf
end