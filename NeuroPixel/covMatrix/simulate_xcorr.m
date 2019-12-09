nsim = 20000; %numbers of simulation
ma=zeros(1,nsim);
maxlag = 10;

n1=10; %number of bins with spikes in v1
n2=1; %number of bins with spikes on v1

for in = 1:nsim
v1=zeros(200,1);
v2=v1;
%number of bins with spikes for vector 1 and 2

%generate a random vector that contains n bins with spikes, each of these
%bins has a random number of spikes
idx = randi(200,[n1,1]);
vals = randi(10,[n1,1]);
%vals = ones([n1,1]);
v1(idx)=vals;
% a second one
idx = randi(200,[n2,1]);
vals = randi(10,[n2,1]);%
%vals = ones(n2,1);
v2(idx)=vals;
%some filtering (not sure that matters)
filt = gausswin(11);
filt = filt/sum(filt);

v1=conv(v1,filt,'same');
v2=conv(v2,filt,'same');

v1=v1-mean(v1);
v2=v2-mean(v2);

[a,b]=xcorr(v1,v2,maxlag,'coeff');
ma(in)=max(a);
end
figure
histogram(ma,100)
%%

nsim = 20000; %numbers of simulation
ma=zeros(1,nsim);
maxlag = 10;

n1=10; %number of bins with spikes in v1
n2=1; %number of bins with spikes on v1

for in = 1:nsim
v1=zeros(200,1);
v2=v1;
%number of bins with spikes for vector 1 and 2
n1=randi(10);
n2=randi(8);
%generate a random vector that contains n bins with spikes, each of these
%bins has a random number of spikes
idx = randi(200,[n1,1]);
vals = randi(10,[n1,1]);
%vals = ones([n1,1]);
v1(idx)=vals;
% a second one
idx = randi(200,[n2,1]);
vals = randi(10,[n2,1]);%
%vals = ones(n2,1);
v2(idx)=vals;
%some filtering (not sure that matters)
filt = gausswin(11);
filt = filt/sum(filt);

v1=conv(v1,filt,'same');
v2=conv(v2,filt,'same');

v1=v1-mean(v1);
v2=v2-mean(v2);

[a,b]=xcorr(v1,v2,maxlag,'coeff');
ma(in)=max(a);
end
figure
histogram(ma,100)
