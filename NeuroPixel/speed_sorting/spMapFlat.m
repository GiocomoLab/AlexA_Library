spMapFlat = zeros(212,200*208);
for iT=1:212
    for iC=1:208
        tmp = squeeze(spMap(:,iT,iC));
        idx = (iC-1)*200+1:iC*200;
        spMapFlat(iT,idx)=tmp;
    end
end

regs = unique(region