X=zeros(numel(post)-1,nnz(sp.cgs==2));
d_pos = discretize(posx,bins);
for iC=1:length(good_cells)
    spike_idx = sp.clu==good_cells(iC);
    c=histcounts(sp.st(spike_idx),post);
    X(:,iC)=c';
end
%%

%Mdl = fitglm(X(1:10000,:),d_pos(1:10000));
%Mdl = mnrfit(X(1:10000,:),d_pos(1:10000));
Mdl=fitcecoc(X(1:10000,:),d_pos(1:10000));
Yhat = predict(Mdl,X);
%%
figure
plot(Yhat)
hold on
plot(d_pos)
yyaxis('right')
plot(trial_contrast(trial))
%%
error = sqrt((Yhat-d_pos(1:end-1)).^2);