files = dir('F:/temp/peak_list/np*.mat');
PEAKS=[];
PEAKValue=[];
QUANT = [];
DEPTH = [];
centers = -125:-250:-3000;
ACG_LIST = {numel(centers,1)};
for iL = 1:numel(centers);
    ACG_LIST{iL}=[];
end
for iF=1:numel(files)
    load(fullfile(files(iF).folder,files(iF).name))
    for iC=1:numel(peak_list)
        if ~isempty(peak_list(iC).peak_loc);
        PEAKS(end+1)=peak_list(iC).peak_loc(1);
        PEAKValue(end+1)=peak_list(iC).peak_val(1);
        QUANT(end+1)=peak_list(iC).quantile(1);
        DEPTH(end+1)=mec_depth(iC);
        [~,icenter] =min(abs(mec_depth(iC)-centers));
        if peak_list(iC).peak_val(1)>peak_list(iC).quantile(1)
        ACG_LIST{icenter}=cat(1,ACG_LIST{icenter},ACG(iC,:));
        end
        end
        
    end
        
end
%%
figure
histogram(PEAKS(PEAKValue>(QUANT+.00)))
%%
figure
idx = PEAKValue>(QUANT+.01);
scatter(PEAKS(idx),DEPTH(idx))
%%
figure
idx = PEAKValue>(QUANT+.01);
h = histogram2(PEAKS(idx)*10,DEPTH(idx),[12 12],'FaceColor','flat');
%%
p=zeros(12,101);
for iC=1:12
    p(iC,:)=mean(ACG_LIST{iC});
end
figure
imagesc(p,[0 0.3])
