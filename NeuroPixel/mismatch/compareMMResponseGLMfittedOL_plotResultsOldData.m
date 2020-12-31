matfiles = dir('F:\Alex\glmFits_oldData2\*.mat');
MM_All =[];
MM_Model_All=[];
TC_All=[];
opt = load_mismatch_opt;
for iF=1:numel(matfiles)
    
    data=load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    %valid_idx = startsWith(data.reg,'MEC');
    valid_idx = true(size(data.mm_resp,1),1);
    if nnz(valid_idx)==0
        continue
    end
    MM = data.mm_resp(valid_idx,:);
    MM_smooth = smoothdata(MM,2,'gaussian',5);
    MM_smooth = MM_smooth-nanmean(MM_smooth(:,opt.time_vecs<-.1 & opt.time_vecs>-1),2);
    count_vecN = MM_smooth;% ./ nanmedian(MM,2);
%count_vecN=data.mm_resp./max(data.mm_resp,[],2);
%model_n = data.mm_predicted./max(data.mm_predicted,[],2);
    model_n = smoothdata(data.mm_predicted(valid_idx,:)/opt.TimeBin);
    model_n=model_n - nanmean(model_n(:,opt.time_vecs<-.1 & opt.time_vecs>-1),2);
opt = load_mismatch_opt;
resp = mean(data.mm_resp(:,opt.time_vecs>0.1 & opt.time_vecs<0.6),2)-mean(data.mm_resp(:,opt.time_vecs<-.2 & opt.time_vecs>-.7),2);
[~,sid]=sort(resp(valid_idx));

model_mat = false(numel(data.glmData),3);

for iC=1:numel(data.glmData)
    if data.glmData(iC).final_pval<0.05
        model_mat(iC,data.glmData(iC).bestModels)=true;
    end
end
model_mat=model_mat(valid_idx,:);
figure('Name',matfiles(iF).name,'Color','w')
subplot(1,3,1)

imagesc(opt.time_vecs,1:numel(sid),count_vecN(sid,:),[-10 10]) 
xlim([-1 3])
title('Actual Response')
subplot(1,3,2)

imagesc(opt.time_vecs,1:numel(sid),model_n(sid,:),[-10 10])
title('Predicted Response')
xlim([-1 3])
subplot(1,3,3)
imagesc(model_mat(sid,:))
set(gca,'XTick',[1:3],'XTickLabel',{'Run','Vis flow','Position'},'XTickLabelRotation',45)
% pause
MM_All = cat(1,MM_All,MM_smooth);
MM_Model_All = cat(1,MM_Model_All,model_n);
tc = nan(nnz(valid_idx),3,200);
tmp = {data.glmData(valid_idx).tuning_curves};
final_pvals = [data.glmData(valid_idx).final_pval];
for iC=1:numel(tmp)
    for iT=1:numel(tmp{iC})
        if final_pvals(iC)<0.05
        if ~isempty(tmp{iC}{iT})
            tc(iC,iT,:)=tmp{iC}{iT};
        end
        end
    end
end
TC_All = cat(1,TC_All,tc);
end

%%
with_tc= ~isnan(TC_All(:,1,1)) & ~isnan(TC_All(:,2,1));
resp = mean(MM_All(with_tc,opt.time_vecs>0.1 & opt.time_vecs<0.6),2);
[~,sid]=sort(resp)
figure('Color','white')
tmp = MM_All(with_tc,:);
tmp_model = MM_Model_All(with_tc,:);
subplot(1,2,1)
imagesc(opt.time_vecs,1:numel(sid),tmp(sid,:),[-10 10])
xlim([-1 3])
title('real response')
subplot(1,2,2)
imagesc(opt.time_vecs,1:numel(sid),tmp_model(sid,:),[-10 10])
xlim([-1 3])
title('predicted response')
%%
figure('Color','White')
TC_tmp = TC_All(with_tc,:,:);
for ii=1:2
    subplot(1,2,ii)
tmp = zscore(squeeze(TC_tmp(sid,ii,:)),1,2);
tmp=tmp(~isnan(tmp(:,1)),:);
imagesc(tmp,[-2 2])
map = brewermap(35,'RdBu');
map = flipud(map);
colormap(map)
end
subplot(1,2,1)
title('Speed Tuning')
subplot(1,2,2)
title('Vis Flow Tuning')
%%

figure
for iC=1:50
    tc=data.glmData(sid(numel(sid)-iC)).tuning_curves;
    
    for iT=1:numel(tc)
        subplot(1,5,iT)
        if ~isempty(tc{iT})
            plot(tc{iT})
        end
    end
    pause
    clf
end