load('Z:\giocomo\attialex\Images\xcorrv_range_gains\allData_MEC_0.6_100.mat')
ff=nanmean(output.XCORRS,3);
figure
hold on
cmap =[0,0,0; cool(4)];

for ii=1:5
plot([-20:2:21],ff(ii,:),'Color',cmap(6-ii,:))
end
xlabel('shift [cm]')
ylabel('XCorr')
legend({'0.5','0.6','0.7','0.8','Baseline'})
grid on

%%
gains = [1.0, 0.8, 0.7,0.6, 0.5];
PEAKS = [];
SHIFTS = [];
cmap =[0,0,0; cool(4)];

for iS=1:23
    figure
    hold on
for iG = 1:5
    tidx = output.GAINS(:,iS)==gains(iG);
    plot(find(tidx),squeeze((output.PEAKS(tidx,:,iS))),'.','Color',cmap(iG,:),'MarkerSize',12);

end

end
    
%%
gains = [1.0, 0.8, 0.7,0.6, 0.5];
PEAKS = [];
SHIFTS = [];
cmap =[0,0,0; cool(4)];

for iS=1:23
    tmpP = zeros(1,5);
    tmpS = zeros(1,5);
for iG = 1:5
    tidx = output.GAINS(:,iS)==gains(iG);
    tmpP(iG) = squeeze(nanmean(output.PEAKS(tidx,:,iS),1));
    tmpS(iG) = squeeze(nanmean(output.SHIFTS(tidx,:,iS)*2,1));
end
PEAKS = cat(1,PEAKS,tmpP);
SHIFTS = cat(1,SHIFTS,tmpS);
end
    
figure
hold on
for ii=1:5
    x=nanmean(SHIFTS(:,ii));
    y=nanmean(PEAKS(:,ii));
    xpos=nanstd(SHIFTS(:,ii))/sqrt(23);
    xneg = xpos;
    ypos = nanstd(PEAKS(:,ii))/sqrt(23);
    yneg = ypos;
    errorbar(x,y,yneg,ypos,xneg,xpos,'o','Color',cmap(ii,:))
end
legend({'Baseline','0.8','0.7','0.6','0.5'})
legend('Location','NW')
xlabel('shift [cm]')
ylabel('XCorr')
grid on
   %%
   figure
   hold on
   for iT = 1:4
   gains = [0.8, 0.7,0.6, 0.5];
   cmap = cool(4);
PEAKS = [];
SHIFTS = [];
for iS=1:23
    tmpP = zeros(1,4);
    tmpS = zeros(1,4);
for iG = 1:4
    tidx = output.GAINS(:,iS)==gains(iG);
    tmpidx = strfind(tidx',[0 1])+iT;
    tmpP(iG) = squeeze(nanmean(output.PEAKS(tmpidx,:,iS),1));
    tmpS(iG) = squeeze(nanmean(output.SHIFTS(tmpidx,:,iS)*2,1));
end
PEAKS = cat(1,PEAKS,tmpP);
SHIFTS = cat(1,SHIFTS,tmpS);
end
    

hold on
for ii=1:4
    x=nanmean(SHIFTS(:,ii));
    y=nanmean(PEAKS(:,ii));
    xpos=nanstd(SHIFTS(:,ii))/sqrt(23);
    xneg = xpos;
    ypos = nanstd(PEAKS(:,ii))/sqrt(23);
    yneg = ypos;
    errorbar(x,y,yneg,ypos,xneg,xpos,'o','Color',cmap(ii,:))
    text(x,y,num2str(iT))
end
   end
   legend({'0.8','0.7','0.6','0.5'})
legend('Location','NW')
legend('Location','NW')
xlabel('shift [cm]')
ylabel('XCorr')
grid on
grid on