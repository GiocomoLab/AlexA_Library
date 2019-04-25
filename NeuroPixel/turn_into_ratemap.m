good_cells = sp.cids(sp.cgs==2);


valid_idx=ismember(sp.clu,good_cells);

st=round(sp.st(valid_idx)*1000)+1;
rastermap=zeros(length(good_cells),max(st));
c_id=sp.clu(valid_idx)+1;

%map cluid to idx
mapping(good_cells+1)=1:length(good_cells);
c_idx=mapping(c_id);
for ii=1:length(c_idx)
    rastermap(c_idx(ii),st(ii))=1;
end
%%
kernel = reshape(gausswin(501),1,[]);
kernel = kernel/sum(kernel);
rate_mat = conv2(rastermap,kernel,'same');
rate_mat = rate_mat/0.001;
rate_mat_ds=rate_mat(:,1:10:end);