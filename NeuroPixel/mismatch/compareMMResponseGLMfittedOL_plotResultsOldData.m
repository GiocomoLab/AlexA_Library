matfiles = dir('F:\Alex\glmFits_oldData\*.mat');
MM_All =[];
MM_Model_All=[];
TC_All=[];
for iF=1:numel(matfiles)
    
    data=load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    valid_idx = startsWith(data.reg,'VISp');
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

model_mat = false(numel(data.glmData),5);

for iC=1:numel(data.glmData)
    if data.glmData(iC).final_pval<0.05
        model_mat(iC,data.glmData(iC).bestModels)=true;
    end
end
model_mat=model_mat(valid_idx,:);
% figure('Name',matfiles(iF).name)
% subplot(1,3,1)
% 
% imagesc(opt.time_vecs,1:numel(sid),count_vecN(sid,:),[-10 10]) 
% 
% subplot(1,3,2)
% imagesc(opt.time_vecs,1:numel(sid),model_n(sid,:),[-10 10])
% subplot(1,3,3)
% imagesc(model_mat(sid,:))
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
resp = mean(MM_All(:,opt.time_vecs>0.1 & opt.time_vecs<0.6),2);
[~,sid]=sort(resp)
figure
imagesc(MM_All(sid,:),[-10 10])
%%
figure;
imagesc(zscore(squeeze(TC_All(sid,1,:)),1,2),[-2 2])
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