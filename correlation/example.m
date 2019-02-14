good_cells = sp.cids(sp.cgs==2);
idx=ismember(sp.clu,good_cells):
[a,b]=CCG(sp.st(idx),double(sp.clu(idx)),'binSize',[0.001],'duration',[0.5]);

mex -O CCGHeart.c
mono=bz_MonoSynConvClick(double(tempSP),sp.st(idx),'plot',false);



for ii=1:38;figure;plot(mono.ccgR(:,mono.sig_con(ii,1),mono.sig_con(ii,2)));end
