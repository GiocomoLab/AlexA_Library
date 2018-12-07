speed = diff(posx);
%figure;plot(posx)
jumps=find(speed<-40);
for ii=1:length(jumps)
    speed(jumps(ii))=.5*speed(jumps(ii)-1)+.5*speed(jumps(ii)+1);
end
[sp_Mat,zn_all,zn_all_d]=get_spike_mat(sp);
%%
aa=corr(zn_all_d',speed(1:end));

%%


[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
depth=zeros(1,size(sp_Mat,1));
for iC=1:size(sp_Mat,1)
    depth(iC)=mean(spikeDepths(sp.clu==iC-1));
end
%%
max_fr=max(zn_all,[],2);
mean_fr=mean(zn_all_d,2);
figure('Name',filenames{iF})
subplot(1,3,1)
plot(aa,depth,'.')
hold on
plot(aa(sp.cgs==2),depth(sp.cgs==2),'r.')
xlabel('correlation with speed')
ylabel('depth on probe')
legend({'MUA','Good'})
xlim([-.5 .5])

subplot(1,3,2)
plot(max_fr,depth,'.')
hold on
plot(max_fr(sp.cgs==2),depth(sp.cgs==2),'r.')
xlabel('max Firing Rate')
ylabel('depth on probe')
legend({'MUA','Good'})

subplot(1,3,3)
plot(mean_fr,depth,'.')
hold on
plot(mean_fr(sp.cgs==2),depth(sp.cgs==2),'r.')
xlabel('mean Firing Rate')
ylabel('depth on probe')
legend({'MUA','Good'})
%%
max_fr=max(zn_all,[],2);
figure('Name',filenames{iF})
plot(aa,max_fr,'.')
hold on
plot(aa(sp.cgs==2),max_fr(sp.cgs==2),'r.')
xlabel('correlation with speed')
ylabel('Max Firing Rate')
legend({'MUA','Good'})
set(gcf,'Position',[680   229   450   749])
%%
aa1=corr(zn_all_d(:,1:82000)',speed(1:82000));
aa2=corr(zn_all_d(:,82001:end)',speed(82001:end));


figure('Name',filenames{iF})
subplot(1,2,1)
plot(aa1,depth,'.')
hold on
plot(aa1(sp.cgs==2),depth(sp.cgs==2),'r.')
xlabel('correlation with speed')
ylabel('depth on probe')
legend({'MUA','Good'})
xlim([-.5 .5])

subplot(1,2,2)
plot(aa2,depth,'.')
hold on
plot(aa2(sp.cgs==2),depth(sp.cgs==2),'r.')
xlabel('correlation with speed')
ylabel('depth on probe')
legend({'MUA','Good'})
xlim([-.5 .5])
