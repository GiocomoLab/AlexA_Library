function [XTX,time_vec]=getCovMatrix2(data,good_cells,trials,time_bin)


time_range=data.post(ismember(data.trial,trials));
start_time = min(time_range);
stop_time = max(time_range);
time_vec = start_time:time_bin:stop_time;
%%
X = zeros(numel(good_cells),numel(time_vec)-1);
%%
for iC=1:numel(good_cells)
    X(iC,:)=histcounts(data.sp.st(data.sp.clu==good_cells(iC)),time_vec);
end
%%
filt=gausswin(10)';
filt = filt/sum(filt);
X=conv2(X,filt,'same');

%%
X=X-mean(X,2);
%X=normc(X);
normc_fcn = @(m) sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
X=normc_fcn(X);
%%
XTX = X'*X;
