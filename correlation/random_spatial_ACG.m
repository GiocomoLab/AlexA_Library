if ispc()
    data_dir = 'Z:\giocomo\attialex\NP_DATA\';
    image_dir = 'Z:\giocomo\attialex\images\';
else
    data_dir = '/oak/stanford/groups/giocomo/attialex/NP_DATA';
    image_dir = '/oak/stanford/groups/giocomo/attialex/images';
end
DATA_RANDOM=struct;
files = dir(fullfile(data_dir,'npI1_0417_dark_1.mat'));

for iF=1:length(files)
    %load(fullfile(data_dir,files(iF).name));
    
    dp2=diff(posx);
    teleport = find(dp2<-1);
    dp2(teleport)=0.5*dp2(teleport-1)+0.5*dp2(teleport+1);
    distance = [0; cumsum(dp2)];
    
    if min(diff(distance))<-1
        %seems some teleports happen across two frames, go again
        dp2=diff(posx);
        teleport = find(dp2<-1);
        double_teleports = teleport((diff(teleport)==1));
        teleport(find(diff(teleport)==1)+1)=[];
        teleport(ismember(teleport,double_teleports))=[];
        dp2(teleport)=0.5*dp2(teleport-1)+0.5*dp2(teleport+1);
        for dp=1:length(double_teleports)
            tmp_i=double_teleports(dp);
            dp2([tmp_i, tmp_i+1])= 0.5*dp2(tmp_i-1)+0.5*dp2(tmp_i+2);
        end
        distance = [0; cumsum(dp2)];
        %single_teleports = teleport(diff(teleport)>1);
    end
    figure;plot(distance)
    
    binedges = 0:2:max(distance)+2;

    time_per_bin = histcounts(distance, binedges);
    ACG_RANDOM={};
    time_per_bin = time_per_bin * params.TimeBin;
    
    good_cells = sp.cids(sp.cgs==2);
    PXX = {};
    for k=1:length(good_cells)
        cluid = good_cells(k);
        spike_id = sp.clu==cluid;
        if nnz(spike_id)>1
        tmp_acg = zeros(500,4001);
        win = gausswin(31);
        gauss_filter = win/sum(win);
        
        [mean_s,q_s]=shuffled_acg(sp,post,distance,spike_id,binedges,gauss_filter,time_per_bin,500);
        ACG_RANDOM{k}=mean_s;
        ACG_QUANTILES{k} = q_s;
        
        %PXX{k}=pxx;
    end
    ACG_RANDOM=cat(1,ACG_RANDOM{:});
    %PXX=cat(2,PXX{:});
    %DATA_RANDOM(iF).PXX=PXX;
    DATA_RANDOM(iF).freq=f;
    DATA_RANDOM(iF).ACG_QUANTILES = ACG_QUANTILES;
    DATA_RANDOM(iF).ACG = ACG_RANDOM;
    DATA_RANDOM(iF).loc = spacing;
%    PXX_tmp = PXX./sum(PXX);
    pp = 1./f;
    %idx = pp>512;
%     idx = (f>1/512 & f<=1/32);
%     PXX_tmp(~idx,:)=0;
%     [aa,ii]=max(PXX_tmp,[],1);
    idx = f>8.5449e-04;
idx(end)=false;
start = find(idx,1);
start=start-1;
locs = zeros(2,size(PXX,2));
for k=1:length(good_cells)
    tmp = PXX(:,k);
    med = median(tmp);
    [a,b]=sort(PXX(idx,k),'descend');
    diff_left = a-PXX(b+start-1,k);
    diff_right = a-PXX(b+start+1,k);
   
    pot_idx = find(diff_left>0 & diff_right > 0,2);
    if numel(pot_idx)>1
    maxloc = b(pot_idx(1))+start;
    sec_loc = b(pot_idx(2))+start;
    locs(:,k)=[maxloc,sec_loc];
    else
        locs(:,k)=1;
    end
end
    pp = 1./f;
    h=figure();
    histogram(pp(locs(1,:)),0:10:1000)
    %saveas(h,fullfile(image_dir,[files(iF).name '_periodHistogram.png']));
    %close(h)
end
%%


idx = (DATA(1).loc*2)>25 & (DATA(1).loc*2)<2000;
start = find(idx,1);
start=start-1;
max_locations = [];
for ii=1:length(good_cells)
    ACG_RANDOM=DATA(1).ACG(ii,:);
    
    [a,b]=sort(ACG_RANDOM(idx),'descend');
    diff_left = a-ACG_RANDOM(b+start-1);
    diff_right = a-ACG_RANDOM(b+start+1);
    
    pot_idx = find(diff_left>0 & diff_right > 0,1);
    has_min = nnz(ACG_RANDOM(idx)<ACG_QUANTILES{ii}(2,idx))>0;
    pot_idx = b(pot_idx(1))+start;
    above_noise = ACG_RANDOM(pot_idx)>ACG_QUANTILES{ii}(1,pot_idx);
    if ~isempty(pot_idx) && has_min && above_noise
        
        max_locations = cat(1,max_locations,pot_idx);
        
    end
    
end