%%

%files = dir('Z:\giocomo\attialex\distance_coding\data\*.mat');
files = dir('/oak/stanford/groups/giocomo/attialex/distance_coding2/data/*.mat');
data_dir = '/oak/stanford/groups/giocomo/attialex/NP_DATA';
%data_dir = 'F:\NP_DATA';
%figure('Position',[680   321   826   657],'visible',false)

p=gcp('nocreate');
if isempty(p)
    parpool(12);
end

f_vec = 0:0.5:500;
parfor iF=1:numel(files)

    xc_data = matfile(fullfile(files(iF).folder,files(iF).name),'Writable',true);
    if ismember('ACG_THETA',fieldnames(xc_data))
        continue
    end
    
    [~,sid]=sort(xc_data.mec_depth,'descend');
    v_idx =1./xc_data.f_vec > 25 &  1./xc_data.f_vec<600;
    data = load(fullfile(data_dir,files(iF).name));

selected_cells = data.sp.cids(data.sp.cgs==2);
st=data.sp.st;

one_vec = ones(size(st));
ACG_TIME=nan(numel(selected_cells),251);

for ii=1:nnz(selected_cells)
    idx = data.sp.clu==selected_cells(ii);
    if nnz(idx)>5
    tmp = CCG(st(idx),one_vec(idx),'binSize',[0.001],'duration',[0.5]);
    tmp = tmp/max(tmp);
    ACG_TIME(ii,:)=squeeze(tmp(251:end));
    end
end


% 
% subplot(1,3,1)
% ACG=xc_data.ACG;
% imagesc(ACG(sid,:),[0 0.3])
% title('Spatial ACG')
% xlabel('Bin (5cm)')
% ylabel('V->D')
% subplot(1,3,2)
% imagesc(ACG_TIME(sid,:),[0 1])
% title('Temporal ACG')
% xlabel('Time [ms]')
% ylabel('V->D')

xc_data.ACG_TIME = ACG_TIME;
xc_data.tvec = 0:0.001:0.25;

selected_cells = data.sp.cids(data.sp.cgs==2);
selected_idx = ismember(data.sp.clu,selected_cells);
[~,~,selected_clu] = unique(data.sp.clu(selected_idx));
spike_times = double(data.sp.st(selected_idx));
time_bins = 0.001;
discrete_time = round(spike_times/time_bins)+1;
spikeMat = zeros(numel(selected_cells),(ceil(max(spike_times))/time_bins)+1,'single');
% populate matrix
for iS=1:numel(spike_times)
    time_idx=discrete_time(iS);
    u_idx = selected_clu(iS);
    spikeMat(u_idx,time_idx)=spikeMat(u_idx,time_idx)+1;

    
end
[PxxSpikes,FSpikes] = pwelch(spikeMat',[],[],f_vec,1/time_bins);
%[PxxSpikes,FSpikes] = pwelch(spikeMat(1:half,:)',[],[],[],1/time_bins);

% get power in theta range, normalized by other frequencies
theta_range=[5 10];
theta_idx = FSpikes>theta_range(1) & FSpikes<=theta_range(2);
rest_idx = ~theta_idx;
thetaPower = mean(PxxSpikes(theta_idx,:));
restPower = mean(PxxSpikes(rest_idx,:));
rest_idx2 = (FSpikes<theta_range(1) | FSpikes > theta_range(2)) &FSpikes<12.5;
restPower2=mean(PxxSpikes(rest_idx2,:));

thetaPowerN = thetaPower./restPower;
thetaPowerN2 = thetaPower./restPower2;
firing_rate = sum(spikeMat,2)/max(spike_times);
xc_data.normalized_theta = thetaPowerN;
xc_data.normalized_theta2 = thetaPowerN2;
xc_data.theta = thetaPower;
xc_data.firing_rate = firing_rate;
% subplot(1,3,3)
% %plot(thetaPowerN,xc_data.mec_depth,'.')
% plot(thetaPowerN(sid),1:numel(sid),'.')
% set(gca, 'YDir','reverse')
% xlim([0 3])
% ylim([1 numel(sid)])






%pause

%saveas(gcf,sprintf('/oak/stanford/groups/giocomo/attialex/distance_coding2/acg_images/%s.png',files(iF).name))





clf
end

