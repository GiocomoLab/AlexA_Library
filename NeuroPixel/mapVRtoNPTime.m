transitions = find(diff([0; bc; 0]));  
%transitions = find(diff([0; bc; 0])>0);  

transitionsVR = find(diff([0; barcode1; 0]));  
%transitionsVR = find(diff([0; barcode1; 0])>0);
%n_UpsVR = transitionsVR(2:2:end) - transitionsVR(1:2:end);


timeVec=nan(size(bc));

for ii=0:40000

%get VR transition n and n-1
trVR1=transitionsVR(end-ii);
trVR2=transitionsVR(end-(ii+1));



%get NP transition n and n-1
trNP1=transitions(end-ii);
trNP2=transitions(end-(ii+1));
%interpolate 
xx=linspace(0,1,trVR1-trVR2+1);
xq=linspace(0,1,trNP1-trNP2+1);
sampleTime=interp1(xx,timestamp(trVR2:trVR1),xq);

timeVec(trNP2:trNP1)=sampleTime;
end
%%
samples=transitions(45000):transitions(45100);
xx=timeVec(samples);

figure
plot(xx,bc(samples));

%find vr frames
gg=find(timestamp>=min(xx) & timestamp<=max(xx));

xxVR=timestamp(gg);
hold on
plot(xxVR,barcode1(gg))

