good_cells = sp.cids(sp.cgs==2);
idx=ismember(sp.clu,good_cells):
[a,b]=CCG(sp.st(idx),double(sp.clu(idx)),'binSize',[0.001],'duration',[0.5]);
