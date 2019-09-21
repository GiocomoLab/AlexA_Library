%load('F:\G5\1207_mismatch_1\1207_mismatch_1.mat')
speed_t=0.05;
% figure('Name',filenames{iF});; plot(speed)
%
if size(mismatch_trigger,1) ~=1
    mismatch_trigger=mismatch_trigger';
end
if nnz(mismatch_trigger==1)>nnz(mismatch_trigger==0)
    %in some old files mismatch_trigger was actually the move variable,
    %i.e. move ==0 is mismatch
    mismatch_trigger = mismatch_trigger<0.1;
end
speed=true_speed';

all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
run_periods=smooth(speed,25)>speed_t;
run_window=-30:30;
possibles=strfind(run_periods',ones(1,length(run_window)))+floor(.5*length(run_window));


mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
tmp = run_periods' & mismatch_trigger<0.9;
possibles_random = strfind(tmp,ones(1,length(run_window)))+floor(.5*length(run_window));
possibles_random = randsample(possibles_random,200);
if all(size(post)==size(speed))
    speed = speed';
end
%% MM
[spike_mat,win,adata]=extract_triggered_spikes(sp,post(mm_trigs),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);
[spike_mat_random,~,adata_random]=extract_triggered_spikes(sp,post(possibles_random),'win',[-4 4],'aux',[post'; [speed]],'aux_win',[-200 200]);

spike_mat=spike_mat(sp.cids+1,:,:);
spike_mat_random=spike_mat_random(sp.cids+1,:,:);
%% get downsampled ratemaps

stepsize = 20;
rate_mat = zeros(size(spike_mat,1),size(spike_mat,2),floor(8001/stepsize));
rate_mat_random = zeros(size(spike_mat_random,1),size(spike_mat_random,2),floor(80001/stepsize));
cntr = 0;
for iS = 1:stepsize:(8001-stepsize)
    cntr = cntr+1;
    idx = iS:iS+stepsize;
    n_spikes =  sum(spike_mat(:,:,idx),3);
    n_spikes_random = sum(spike_mat_random(:,:,idx),3);
    rate_mat(:,:,cntr)=n_spikes;
    rate_mat_random(:,:,cntr)=n_spikes_random;
end
%%
    tmpMM=squeeze(mean(rate_mat,2));
    mmresp = mean(tmpMM(:,220:250),2)-mean(tmpMM(:,155:185),2);
    
    [~,sidx]=sort(mmresp,'descend');

%% calculate significant bins
sig_val = zeros(size(rate_mat,1),size(rate_mat,3));


for iC = 1:size(rate_mat,1)
for iS=1:size(sig_val,2)



    [h,p]=ttest2(squeeze(rate_mat_random(iC,:,iS)),squeeze(rate_mat(iC,:,iS)));
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
frac = .05;
n=round(frac*size(rate_mat,1));
figure
tvec = linspace(-4,4,size(rate_mat,3));
tmp = mean(squeeze(mean(rate_mat(sidx(1:n),:,:),2)));
plot(tvec,tmp)
hold on
ff=nanmedian(crossings(sidx(1:n)));
ff=round(ff)+win(1)-1;
plot(tvec(ff),tmp(ff),'r*')

tmp_c = crossings(sidx(1:n));
tmp_c(isnan(tmp_c))=[];
tmp_c=tmp_c+win(1)-1;
plot(tvec(tmp_c),tmp(ff)*ones(1,numel(tmp_c)),'g.')