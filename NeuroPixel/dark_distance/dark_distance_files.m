if ispc()
    data_dir = 'Z:\giocomo\attialex\NP_DATA\';
    image_dir = 'Z:\giocomo\attialex\images\';
else
    data_dir = '/oak/stanford/groups/giocomo/attialex/NP_DATA';
    image_dir = '/oak/stanford/groups/giocomo/attialex/images';
end
files = dir(fullfile(data_dir,'*dark*.mat'));
params = readtable('UniversalParams.xlsx');
for iF=1:length(files)
    load(fullfile(data_dir,files(iF).name));
    DATA=struct;
    dp2=diff(posx);
    teleport = find(dp2<-1);
    dp2(teleport)=0.5*dp2(teleport-1)+0.5*dp2(teleport+1);
    total_distance = [0; cumsum(dp2)];
    
    if min(diff(total_distance))<-1
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
        total_distance = [0; cumsum(dp2)];
        %single_teleports = teleport(diff(teleport)>1);
    end
    figure;plot(total_distance)
    
    binedges = 0:2:max(total_distance)+2;
    
    time_per_bin = histcounts(total_distance, binedges);
    ACG={};
    ACG_RANDOM={};
    time_per_bin = time_per_bin * params.TimeBin;
    
    good_cells = sp.cids(sp.cgs==2);
    PXX = {};
    for k=1:length(good_cells)
        cluid = good_cells(k);
        spike_id = sp.clu==cluid;
        spike_t = sp.st(spike_id);
        [~,~,spike_idx] = histcounts(spike_t,post);
        spike_total_distance = total_distance(spike_idx);
        [aa,bins]=histcounts(spike_total_distance,binedges);
        
        %moothing
        firing_rate = aa./time_per_bin;
        if sum(isnan(firing_rate))>0
            firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
        end
        % smoothSigma = 3;%params.SmoothSigmaFR/params.SpatialBin;
        % smoothWindow = floor(smoothSigma*5/2)*2+1;
        % gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
        win = gausswin(31);
        gauss_filter = win/sum(win);
        %aa = conv(firing_rate,gauss_filter,'same');
        
        [acg,spacing,pxx,f] = calculate_spatial_acg(spike_total_distance,binedges,gauss_filter,time_per_bin);
        [mean_s,q_s]=shuffled_acg(sp,post,total_distance,spike_id,binedges,gauss_filter,time_per_bin,500);
        ACG_RANDOM{k}=mean_s;
        ACG_QUANTILES{k} = q_s;
        %[acg,spacing] = xcorr(aa,2000,'coeff');
        %[pxx,f]=pwelch(aa,[],[],[],1/2);
        ACG{k}=acg;
        PXX{k}=pxx;
        %[shuffled_acg,~]=
    end
    ACG=cat(1,ACG{:});
    PXX=cat(2,PXX{:});
    ACG_RANDOM = cat(1,ACG_RANDOM{:});
    DATA.PXX=PXX;
    DATA.freq=f;
    DATA.ACG = ACG;
    DATA.ACG_RANDOM = ACG_RANDOM;
    DATA.ACG_QUANTILES = ACG_QUANTILES;
    DATA.loc = spacing;

    save(fullfile(data_dir,files(iF).name),'DATA','total_distance','-append');
    
end

%%

all_freqs=[];
for iF=1:length(DATA)
    idx = DATA(iF).freq>8.5449e-04;
    PXX = DATA(iF).PXX;
    idx(end)=false;
    start = find(idx,1);
    start=start-1;
    locs = zeros(2,size(PXX,2));
    for k=1:size(PXX,2)
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
    all_freqs=cat(1,all_freqs,DATA(iF).freq(locs(2,:)));
end

histogram(all_freqs,0:0.001:0.03)

%% same for ccg
all_Locations=[];
for iF=1:length(DATA)
    idx = (DATA(iF).loc*2)>25 & (DATA(iF).loc*2)<1000;
    start = find(idx,1);
    start=start-1;
    ACG=DATA(iF).ACG;
    locs = zeros(2,size(ACG,1));
    
    for k=1:size(ACG,1)
        tmp = ACG(k,:);
        med = median(tmp);
        [a,b]=sort(ACG(k,idx),'descend');
        diff_left = a-ACG(k,b+start-1);
        diff_right = a-ACG(k,b+start+1);
        
        pot_idx = find(diff_left>0 & diff_right > 0,2);
        
        if numel(pot_idx)>1
            maxloc = b(pot_idx(1))+start;
            sec_loc = b(pot_idx(2))+start;
            locs(:,k)=[maxloc,sec_loc];
        else
            locs(:,k)=1;
        end
    end
    
    
    all_Locations=cat(2,all_Locations,DATA(iF).loc(locs(1,:)));
end
figure
histogram(all_Locations,200)
%%
idx = (DATA(1).loc*2)>25 & (DATA(1).loc*2)<2000;
start = find(idx,1);
start=start-1;
max_locations = [];
spacing = DATA(1).loc;
for ii=1:length(good_cells)
    ACG=DATA(1).ACG(ii,:);
    
    [a,b]=sort(ACG(idx),'descend');
    diff_left = a-ACG(b+start-1);
    diff_right = a-ACG(b+start+1);
    
    pot_idx = find(diff_left>0 & diff_right > 0,1);
    has_min = nnz(ACG(idx)<ACG_QUANTILES{ii}(2,idx))>0;
    pot_idx = b(pot_idx(1))+start;
    above_noise = ACG(pot_idx)>ACG_QUANTILES{ii}(1,pot_idx);
    if ~isempty(pot_idx) && has_min && above_noise
        
        max_locations = cat(1,max_locations,spacing(pot_idx)*2);
        
    end
    
end
