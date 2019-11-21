function [ACG_corr,ZERO_LAG_PRE,ZERO_LAG_POST] = extractXCorr(data,trials,region)

%cellIDX=find(sp.cgs>=1);
if isfield(data.anatomy,'parent_shifted')
    reg = startsWith(data.anatomy.parent_shifted,region);
else
reg = startsWith(data.anatomy.cluster_parent,region);
end
if iscolumn(reg)
    reg = reg';
end

selected_cells = data.sp.cids(data.sp.cgs==2 & reg);
st=data.sp.st(ismember(data.sp.clu,selected_cells));
clu = data.sp.clu(ismember(data.sp.clu,selected_cells));
[a,~,clu]=unique(clu);
%
time_range=data.post(ismember(data.trial,trials));
start_time = min(time_range);
stop_time = max(time_range);
duration = stop_time-start_time;

pre_idx = st>start_time-duration & st<start_time;
post_idx = st>start_time & st<stop_time;
[CGR_pre,b]=CCG(st(pre_idx),clu(pre_idx),'binSize',[0.001],'duration',[0.5]);
[CGR_post,b]=CCG(st(post_idx),clu(post_idx),'binSize',[0.001],'duration',[0.5]);

%extract autocorrelation
ACG_PRE = zeros(size(CGR_pre,1),size(CGR_pre,2));
for iC=1:size(CGR_pre,2)
    ACG_PRE(:,iC)=squeeze(CGR_pre(:,iC,iC));
end
ACG_POST = zeros(size(CGR_post,1),size(CGR_post,2));
for iC=1:size(CGR_pre,2)
    ACG_POST(:,iC)=squeeze(CGR_post(:,iC,iC));
end
ACG_corr = diag(corr(ACG_PRE,ACG_POST));
%extract zero lag xcorr for each pair
npairs = size(CGR_pre,2);
npairs = npairs*(npairs-1)/2;
ZERO_LAG_PRE = zeros(npairs,1);
ZERO_LAG_POST = ZERO_LAG_PRE;
midpoint = 251;
cntr = 0;
n_cells = numel(a);
for ii =1:numel(a)
    for jj = (ii+1):numel(a)
        cntr=cntr+1;
        %norm_pre = squeeze(sqrt(CGR_pre(midpoint,ii,ii)*CGR_pre(midpoint,jj,jj)));
        ZERO_LAG_PRE(cntr)=squeeze(CGR_pre(midpoint,ii,jj));
        ZERO_LAG_POST(cntr)=squeeze(CGR_post(midpoint,ii,jj));
    end
end
end

