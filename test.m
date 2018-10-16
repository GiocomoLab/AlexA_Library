PS=dataArray{:,3}~=1;

velM=dataArray{:,2};
figure
plot(velM)
velMBi=velM>0.05;
hold on
plot(velMBi)
%%
trigs=strfind(velMBi',[zeros(1,30) ones(1,30)])+29;

figure
plot(velM)
hold on
plot(trigs,velM(trigs),'ro')

%%
 trigs=strfind(PS',[0 1]);
 running_period=-30:30;
  vidx=false(size(trigs));
                for iT=1:length(trigs)
                    if sum(velMBi(trigs(iT)+running_period))==length(running_period)
                        vidx(iT)=true;
                    end
                end
 
 %%
