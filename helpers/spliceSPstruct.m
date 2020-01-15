function sp_new = spliceSPstruct(sp,tstart,tstop)

n_chunks = numel(tstart);
st_new = [];
clu_new=[];
for ii=1:n_chunks
    idx = sp.st>=tstart(ii) & sp.st<tstop(ii);
    if ii>1
    %tmp = sp.st(idx(1))-tstart(ii);
    tmp = tstart(ii)-st_new(end);
    else
        tmp = tstart(ii);
    end
    st_new = cat(1,st_new,sp.st(idx)-tmp);
    clu_new = cat(1,clu_new,sp.clu(idx));
end

sp_new.st = st_new;
sp_new.clu = clu_new;
end