%matfiles = dir(fullfile(OAK,'attialex','gain_changeMapsLargeBin','*.mat'));
matfiles = dir(fullfile(OAK,'attialex','gain_changeMaps','*.mat'));

decode_error=[];
distance =[];

for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    decode_error=cat(1,decode_error,data_out.score_mat(1,:));
    distance = cat(1,distance,data_out.score_mat(2,:));
end

%%
decode_error(abs(decode_error)>50)=nan;
error_pos = nanmean(decode_error(:,251:280),2);

slices = 4;
CT=cbrewer('qual', 'Set2',slices);
N=floor(numel(error_pos)/slices);
[a,b]=sort(error_pos);
figure
for iS=1:slices
idx = (iS-1)*N+1:iS*N;
idx=b(idx);
tmp_pos = nanmean(decode_error(idx,:));
tmp_distance = nanmean(distance(idx,:));
subplot(2,1,1)
hold on
%boundedline(x_vec,tmp_pos,std(CV(idx,:))/sqrt(nnz(idx)),'cmap',CT(iS+2,:))

plot(tmp_pos,'Color',CT(iS,:))
subplot(2,1,2)
hold on
%boundedline(x_vec,tmp_distance,std(SV(idx,:))/sqrt(nnz(idx)),'cmap',CT(iS+2,:))

plot(tmp_distance,'Color',CT(iS,:))
end