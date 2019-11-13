%%
[XTX]=getCovMatrix2(data,good_cells,[11:34],0.25);
%%

figure
imagesc(XTX,[0 1])
hold on
%xline(gain_idx,'r')
%yline(gain_idx,'r')


%%
tt=corr(X(:,1:400)');

ff=corr(X(:,4001:end)');
figure;subplot(1,2,1);imagesc(tt,[-.5 .5]);subplot(1,2,2);imagesc(ff,[-.5 .5])

%%
time_bin = 0.1;
time_range=data.post(ismember(data.trial,[16:20]));
start_time = min(time_range);
stop_time = max(time_range);
time_vec = start_time:time_bin:stop_time;

X = zeros(numel(good_cells),numel(time_vec)-1);

for iC=1:numel(good_cells)
    X(iC,:)=histcounts(data.sp.st(data.sp.clu==good_cells(iC)),time_vec);
end

filt=gausswin(1)';
filt = filt/sum(filt);
X=conv2(X,filt,'same');
X=X-mean(X,2);
%X=normc(X);
normc_fcn = @(m) sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
X=normc_fcn(X);
%%
figure
imagesc(corr(X'))
%%
Y = pdist(X);
Z = linkage(Y);
order = optimalleaforder(Z,Y);
figure
tmp = corr(X');
imagesc(tmp(order,order))


%%

%%
time_bin = 0.1;
time_range=data.post(ismember(data.trial,[21:24]));
start_time = min(time_range);
stop_time = max(time_range);
time_vec = start_time:time_bin:stop_time;

Y = zeros(numel(good_cells),numel(time_vec)-1);

for iC=1:numel(good_cells)
    Y(iC,:)=histcounts(data.sp.st(data.sp.clu==good_cells(iC)),time_vec);
end

filt=gausswin(1)';
filt = filt/sum(filt);
Y=conv2(Y,filt,'same');
Y=Y-mean(Y,2);
%X=normc(X);
normc_fcn = @(m) sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
Y=normc_fcn(Y);

%%
figure
tmp = corr(X');
tmp2 = corr(Y');
subplot(1,2,1)
imagesc(tmp(order,order))
subplot(1,2,2)
imagesc(tmp2(order,order))