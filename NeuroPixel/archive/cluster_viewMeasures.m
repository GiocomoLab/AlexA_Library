%start phy in debugging mode, get text of table:
%txt=controller.supervisor.cluster_view.html()
%write it to html file fn=open('cluter.html','w'), fn.write(txt),
%fn.close()
%open in explorer
%copy paste to excel...
%%
goodIDX=Classified=='good';

figure
histogram(VR(goodIDX),50,'DisplayStyle','stairs','Normalization','probability','BinLimits',[0 0.02])
hold on
histogram(VR(~goodIDX),50,'DisplayStyle','stairs','Normalization','probability','BinLimits',[0 0.02])

xlabel('Fraction of spikes within 2ms')

ylabel('Count')

legend({'Single Units','Rest'})

%%
goodIDX=Classified=='good';

figure
histogram(R(goodIDX),50,'DisplayStyle','stairs','Normalization','probability','BinLimits',[0 1])
hold on
histogram(R(~goodIDX),50,'DisplayStyle','stairs','Normalization','probability','BinLimits',[0 1])

xlabel('ratio 1ms bin to rest')

ylabel('Count')

legend({'Single Units','Rest'})

%%
goodIDX=Classified=='good';

figure
histogram(amplitude(goodIDX),50,'DisplayStyle','stairs','Normalization','probability','BinLimits',[0 150])
hold on
histogram(amplitude(~goodIDX),50,'DisplayStyle','stairs','Normalization','probability','BinLimits',[0 150])

xlabel('Amplitude')

ylabel('Count')

legend({'Single Units','Rest'})
