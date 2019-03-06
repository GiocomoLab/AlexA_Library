function firing_rate = get_firing_rate_glm(A,bestModels,parameters)
%calculates firing rate based on 
X  = ones(length(A{1}),1);
for ib=bestModels
    X=[X, A{ib}];
end
firing_rate = exp(X*parameters);
end
% figure
% hold on
% plot(sp)
% 
% win=gausswin(41);
% win=win/sum(win);
% fr=conv(spiketrain,win,'same');
% plot(fr)