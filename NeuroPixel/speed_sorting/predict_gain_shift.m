ALL_M=[];
ALL_S=[];
ALL_G=[];
ALL_SG=[];
ALL_PRED=[];
REGION = [];
STABILITY=[];
for iS = 1:numel(output)
data_out = output{iS};
nC=numel(data_out.region);
delta_m=zeros(1,nC);
delta_s = delta_m;
delta_g = delta_m;
pred_m = delta_m;
delta_sg=delta_m;
for iC=1:nC
    [~,max_slow]=max(data_out.all_slow(iC,:));
    [~,max_fast] = max(data_out.all_fast(iC,:));
    [~,max_gain] = max(data_out.all_gain(iC,:));
    
    trial_speed_all = getSpeedAroundPoint(data_out.speed(:,1),data_out.speed(:,2),data_out.speed(:,3),ops,data_out.max_ind(iC),ops.speedWindow);
    
    trial_speed = trial_speed_all(1:end-4);
    [speed_sorted,sidx]=sort(trial_speed,'descend');
    
    fast_speed = mean(speed_sorted(1:2));
    slow_speed = mean(speed_sorted(end-1:end));
    gain_speed = mean(trial_speed_all(end-3:end));
    
    delta_m(iC) = max_fast-max_slow;
    delta_s(iC) = fast_speed-slow_speed;
    delta_g(iC)=slow_speed - gain_speed;
    delta_sg(iC) = max_slow-max_gain;
    
    pred_m(iC)= -(fast_speed - slow_speed)*data_out.factors(iC);
end
ALL_M = cat(2,ALL_M,delta_m);
STABILITY = cat(1,STABILITY,data_out.stability);
REGION = cat(1,REGION,data_out.region');
ALL_PRED = cat(2,ALL_PRED,pred_m);
ALL_S = cat(2,ALL_S,delta_s);
ALL_G = cat(2,ALL_G,delta_g);
ALL_SG = cat(2,ALL_SG,delta_sg);
end
%%
figure
idx = STABILITY>.2 & ismember(REGION,'MEC');
plot(ALL_M(idx),ALL_PRED(idx),'.')
xlabel('actual difference between peak loc')
ylabel('predicted difference between peak loc')
%% now do something similar 
figure
plot(ALL_M(idx),ALL_S(idx),'.')
hold on
plot(ALL_SG(idx),ALL_G(idx),'r.')