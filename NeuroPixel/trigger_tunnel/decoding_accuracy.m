%data=load('/Volumes/T7/AA_201025_3_TowerTraining_201125_18-27-26.mat');
accuracy_all = {};
mf = dir('F:\Alex\new_2\*Tower*.mat');
for iF=1:numel(mf)
    data = load(fullfile(mf(iF).folder,mf(iF).name));
good_cells = data.sp.cids(data.sp.cgs==2);
tmp = diff(data.posx)<-50;
data.trial = [0;cumsum(tmp)];
opt = load_mismatch_opt;

%%
trigs_1=strfind(data.vr_data_resampled.Object==1,[0 0 1 1])+2;
trigs_2=strfind(data.vr_data_resampled.Object==0,[0 0 1 1])+2;
[spikeTimes1,~,~,~,count_vec1]=extract_triggered_spikeTimes(data.sp,data.post(trigs_1),'cluIDs',good_cells,'win',opt.extract_win);

[spikeTimes2,~,~,~,count_vec2]=extract_triggered_spikeTimes(data.sp,data.post(trigs_2),'cluIDs',good_cells,'win',opt.extract_win);

%%
stim_1_ind=zeros(numel(good_cells),numel(opt.time_vecs),numel(spikeTimes1));
for ii=1:numel(spikeTimes1)
    stim_1_ind(:,:,ii)= spikeStruct2countVec(spikeTimes1{ii},good_cells,opt);
end

stim_2_ind=zeros(numel(good_cells),numel(opt.time_vecs),numel(spikeTimes2));
for ii=1:numel(spikeTimes2)
    stim_2_ind(:,:,ii)= spikeStruct2countVec(spikeTimes2{ii},good_cells,opt);
end

X=cat(3,stim_1_ind,stim_2_ind);
labels=[ones(numel(spikeTimes1),1) ;zeros(numel(spikeTimes2),1)];

%%
X_new = smoothdata(X,2,'movmean',11);
sample_idx = 1:11:numel(opt.time_vecs);
%X_new = X_new(:,1:5:end,:);
%%
accuracy = zeros(1,size(X_new,2));
loss = accuracy;
for ii=1:numel(accuracy)
    obs = squeeze(X_new(:,ii,:));
    CVMdl = fitclinear(obs,labels,'ObservationsIn','columns','KFold',5,'Learner','logistic','Solver','sparsa','Regularization','lasso');
    yhat = CVMdl.kfoldPredict;
    accuracy(ii) = nnz(yhat==labels)/numel(labels);
    loss(ii)=kfoldLoss(CVMdl);
end
accuracy_all{end+1}=accuracy;
figure('Name',mf(iF).name)
plot(opt.time_vecs,accuracy)

%% refit best model
[~,mi]=max(accuracy);
obs = squeeze(X_new(:,mi,:));
    %[Mdl,fitinfo] = fitclinear(obs,labels,'ObservationsIn','columns','Learner','logistic','Solver','sparsa','Regularization','lasso');
CVMdl = fitclinear(obs,labels,'ObservationsIn','columns','KFold',5,...
    'Learner','logistic','Solver','sparsa','Regularization','lasso');
    fprintf('best model at %.2f, accuracy: %.2f \n',opt.time_vecs(mi),1-kfoldLoss(CVMdl))
    

    %% important cells:
        bias=[];
    for ii=1:numel(CVMdl.Trained);
        bias = cat(2,bias,CVMdl.Trained{ii}.Beta);
    end
    sig_cells = all(abs(bias)>0,2);
    [sig_idx]=find(sig_cells);
    figure('Name',mf(iF).name)
    cv1=smoothdata(count_vec1,2,'gaussian',11);
    cv2=smoothdata(count_vec2,2,'gaussian',11);
    for iC=1:min(numel(sig_idx),12)
        subplot(4,3,iC)
        plot(opt.time_vecs,count_vec1(sig_idx(iC),:))
        hold on
        plot(opt.time_vecs,count_vec2(sig_idx(iC),:))
        

    end
end