%% run loop_mm_across_region first
%% speed tuning
nSlices = 7;
cmap = cbrewer('div','RdYlBu',nSlices);
cmap = flipud(cmap);
cmap(3,:)=[.2 1 .20];
avg_all = zeros(numel(RUN_TRACES),nSlices,numel(opt.time_vecs));
RUN_TRACES = RUN_TRACES(~cellfun('isempty',RUN_TRACES));
CLUIDS = CLUIDS(~cellfun('isempty',CLUIDS));
SPIKE_TIMES = SPIKE_TIMES(~cellfun('isempty',SPIKE_TIMES));
figure('Color','White')
usites = unique(SID);
PARAMS = [];
GOF = [];
REG = {};
for iSite = 1:numel(RUN_TRACES)

avgMM=mean(MM(SID==usites(iSite),opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(SID==usites(iSite),opt.time_vecs<-.1 & opt.time_vecs>-0.5),2);
[a,b]=sort(avgMM,'descend');
mm_increase = CLUIDS{iSite}(b(1:round(numel(b)/5)));

run_speed = mean(RUN_TRACES{iSite}(:,25:50),2);
[avg_speed,run_si]=sort(run_speed);
trial_vec =cat(1,SPIKE_TIMES{iSite}{:});
nTrials = size(RUN_TRACES{iSite},1);
count_vec = zeros(nTrials,numel(opt.time_bins)-1);

for iT=1:nTrials
    idx = trial_vec(:,3)==iT & ismember(trial_vec(:,2),mm_increase);
    [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
    count_vec(iT,:)=spike_count;
end

chunkSize = round(nTrials/nSlices);

hold on
speed=zeros(1,nSlices);
resp = zeros(1,nSlices);
resp_1 = mean(count_vec(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)- mean(count_vec(:,opt.time_vecs>-.5 & opt.time_vecs<-.1),2);

for iS=1:nSlices
    subplot(2,1,1)
    hold on
    idx = (iS-1)*chunkSize +1 : (min(iS*chunkSize,nTrials));
    avg = mean(count_vec(run_si(idx),:)-mean(count_vec(run_si(idx),opt.time_vecs<-.1 & opt.time_vecs>-.5),2));
    avg_all(iSite,iS,:)=avg;
    plot(opt.time_vecs,avg,'Color',cmap(iS,:),'LineWidth',2)
%     subplot(1,2,2)
%     hold on
    %plot(mean(RUN_TRACES{iSite}(run_si(idx),:)))
    speed(iS)=mean(avg_speed(idx));
    resp(iS)=mean(resp_1(run_si(idx)));
    
end
xlim([-.5 1])
xlabel('time')
ylabel('Firing Rate')
subplot(2,1,2)
%speed= avg_speed;
%resp = mean(count_vec(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)- mean(count_vec(:,opt.time_vecs>-.5 & opt.time_vecs<-.1),2);
%scatter(speed,resp,'ro');
hold on
speed_1= run_speed;
%resp_1 = mean(count_vec(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)- mean(count_vec(:,opt.time_vecs>-.5 & opt.time_vecs<-.1),2);
scatter(speed_1,resp_1,'.')
[p,S]=polyfit(speed_1,resp_1,1);
yfit = polyval(p,linspace(min(speed_1),max(speed_1),5));
hold on
plot(linspace(min(speed),max(speed),5),yfit)
reg_this = REGION(SID==usites(iSite));

title(reg_this{1})
%xlim([-.5 1])
xlabel('speed')
ylabel('Firing Rate')
%pause
clf

PARAMS = cat(1,PARAMS,p);
GOF = cat(1,GOF,S.normr);
REG = cat(1,REG,reg_this(1));

end
%%

nSlices = 5;
avg_all = zeros(numel(RUN_TRACES),nSlices,numel(opt.time_vecs));
RUN_TRACES = RUN_TRACES(~cellfun('isempty',RUN_TRACES));
CLUIDS = CLUIDS(~cellfun('isempty',CLUIDS));
SPIKE_TIMES = SPIKE_TIMES(~cellfun('isempty',SPIKE_TIMES));

usites = unique(SID);
PARAMS = [];
GOF = [];
REG = {};
for iSite = 1:numel(RUN_TRACES)

avgMM=mean(MM(SID==usites(iSite),opt.time_vecs>0.1 & opt.time_vecs<0.5),2)-mean(MM(SID==usites(iSite),opt.time_vecs<-.1 & opt.time_vecs>-0.5),2);
reg_this = REGION(SID==usites(iSite));
[a,b]=sort(avgMM,'descend');
mm_increase = CLUIDS{iSite}(b(1:round(numel(b)/5)));

run_speed = mean(RUN_TRACES{iSite}(:,25:50),2);
[avg_speed,run_si]=sort(run_speed);
trial_vec =cat(1,SPIKE_TIMES{iSite}{:});
nTrials = size(RUN_TRACES{iSite},1);
count_vec = zeros(nTrials,numel(opt.time_bins)-1);

for iT=1:nTrials
    idx = trial_vec(:,3)==iT & ismember(trial_vec(:,2),mm_increase);
    [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
    count_vec(iT,:)=spike_count;
end

chunkSize = round(nTrials/nSlices);

speed= run_speed;
resp = mean(count_vec(:,opt.time_vecs>0.1 & opt.time_vecs<0.5),2)- mean(count_vec(:,opt.time_vecs>-.5 & opt.time_vecs<-.1),2);
% speed=zeros(1,nSlices);
% resp = zeros(1,nSlices);
% for iS=1:nSlices
% 
%     idx = (iS-1)*chunkSize +1 : (min(iS*chunkSize,nTrials));
%     avg = mean(count_vec(run_si(idx),:)-mean(count_vec(run_si(idx),opt.time_vecs<0 & opt.time_vecs>-.5),2));
%     avg_all(iSite,iS,:)=avg;
%     speed(iS)=mean(avg_speed(idx));
%     resp(iS)=mean(avg(opt.time_vecs>0.1 & opt.time_vecs<0.6));
% end
[p,S]=polyfit(speed,resp,1);
PARAMS = cat(1,PARAMS,p);
GOF = cat(1,GOF,S.normr);
REG = cat(1,REG,reg_this(1));
end
%%
regions = {'VISp','MEC'};
figure('Color','white')
slope = {};
gof = {};
for iR=1:2
    subplot(1,2,1)
    idx = startsWith(REG,regions{iR});
    slope{iR}=PARAMS(idx,1);
    gof{iR}=GOF(idx);
end
subplot(1,2,1)
plotSpread(slope)
xticklabels(regions)
ylabel('slope')
subplot(1,2,2)
plotSpread(gof)
xticklabels(regions)
ylabel('residuals')