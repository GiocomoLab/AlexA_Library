gains=trial_gain<=.8;
transitions=strfind(gains',[0 1]);

trials=[transitions transitions+1];
spC=[];
dwell_time=[];
for iT=1:length(trials)
    idxVR=trial==trials(iT);
    t_time=post(idxVR);
    start=min(t_time);
    stop=max(t_time);
    idxNP=sp.st<stop & sp.st>=start;
    edges=[0:2:402];
    edges(1)=-5;
    [spM, dT]=getSpikeMatPosition(sp.st(idxNP),sp.clu(idxNP),posx(idxVR),post(idxVR),'edges',edges);
    spC=cat(3,spC,spM);
    dwell_time=cat(1,dwell_time,dT);
end
%cellIDX=find(sp.cgs>=1);
spC=spC(sp.cids+1,:,:);
%%

ffpre=mean(spC(:,:,1:length(transitions)),3);
ffpost = mean(spC(:,:,9:end),3);
mm_pre=mean(ffpre(:,1:5),2);
mm_post=mean(ffpost(:,1:5),2);
figure
plot(mm_pre,mm_post,'.')
grid on
axis image
%%
gains=trial_gain<=.8;
transitions=strfind(gains',[0 1]);

trialsPre=transitions;
trialsPost = transitions+1;
trigs_pre = [];
trigs_post = [];
for ii=1:length(transitions)
    idxPre=find(trial==transitions(ii),1);
    idxPost = find(trial==transitions(ii)+1,1);
    trigs_pre(end+1)=idxPre;
    trigs_post(end+1)=idxPost;
end
win=[-4 4];
[spike_matPre,win,adata,meanPre]=extract_triggered_spikes(sp,post(trigs_pre),'win',win,'aux',[post'],'aux_win',[-200 200]);
[spike_matPost,win,adata,meanPost]=extract_triggered_spikes(sp,post(trigs_post),'win',win,'aux',[post'],'aux_win',[-200 200]);

% binned_array=squeeze(spike_mat(1,:,:));
% bins=win(1):0.001:win(2);
%     
%     % set the new rasters
%     [tr,b] = find(binned_array);
%     [rasterX,yy] = rasterize(bins(b));
%     rasterY = yy*1+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
% figure
% plot(rasterX,rasterY)
%%
meanPre=meanPre(sp.cids+1,:);
meanPost = meanPost(sp.cids+1,:);
figure
plot(mean(meanPre))
hold on
plot(mean(meanPost))
%%
IDX=aggregateData.AID==4;
mmresp=mean(aggregateData.avgMM(IDX,205:250),2)-mean(aggregateData.avgMM(IDX,140:190),2);
gainresp=mean(meanPost(:,4500:6500),2)-mean(meanPre(:,4500:6500),2);
figure
plot(mmresp,gainresp,'.')