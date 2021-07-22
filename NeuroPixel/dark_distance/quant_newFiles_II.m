raw_data_dir = 'F:\Alex\matfiles_new';


%%
visp_files={    
   {'AA_200919_1_MMdark_201010_13-32-24.mat'},{80:160};
    {'AA_200919_1_MMdark_201011_16-11-31.mat'},{85:200};
    {'AA_200919_3_MMdark_201019_12-00-55.mat'},{70:150};
    {'AA_200920_5_MMdark_201015_16-59-26.mat'},{85:180};
    };
 %%
 rs_files = {
    {'AA_200919_1_MMdark_201013_12-03-20.mat'},{130:250};
    {'AA_200919_1_MMdark_201014_12-29-09.mat'},{150:270};
    {'AA_200919_3_MMdark_201015_13-47-57.mat'},{150:300};
    {'AA_200920_5_MMdark_201018_11-48-04.mat'},{150:300};
    {'AA_200920_5_MMdark_201019_14-21-28.mat'},{150:280};
    };
%%
files = cat(1,visp_files,rs_files);
%%
savepath = fullfile('F:\Alex','distance_tuning');

if ~isfolder(savepath)
    mkdir(savepath)
end

%ops = load_default_opt;
ops = load_mismatch_opt;
ops.num_shuf = 400;
ops.dark = true;
ops.SpatialBin = 5;
ops.max_lag = 800;

for iF=1:size(files,1)
    data = load(fullfile(raw_data_dir,files{iF}{1}));
    [~,sn]=fileparts(files{iF}{1});
    clu_id = data.sp.waveform_metrics.cluster_id;
    gc = data.sp.cids(data.sp.cgs==2);
    nC=nnz(ismember(data.sp.waveform_metrics.peak_channel,files{iF,2}{1}) & ismember(clu_id,gc));
    disp(nC)
    data_out = calc_distance_tuning(data,data.sp.cids(data.sp.cgs==2),ops);
    
    save_name = fullfile(savepath,files{iF}{1});

    save(save_name,'data_out')
end


%%
figure
ACG_ALL=[];
DEPTH = [];
ma = 0;
for iF=1:size(visp_files,1)
    data = load(fullfile(savepath,visp_files{iF,1}{1}));
    
    data = data.data_out;
    raw_data = load(fullfile(raw_data_dir,visp_files{iF}{1}));
    
    wft =raw_data.sp.waveform_metrics;
    valid_rows = ismember(wft.cluster_id,data.good_cells);
    chan_number = wft.peak_channel(valid_rows);
    
    idx = ismember(chan_number,visp_files{iF,2}{1});
    tmp = chan_number(idx);
    tmp = tmp-min(tmp);
    DEPTH = cat(1,DEPTH,tmp);
    ACG = data.xcorrs(idx,:);
    %ACG = ACG(sid,:);
    ACG_ALL=cat(1,ACG_ALL,ACG);
end
ma = max(size(ACG_ALL,1),ma);
subplot(1,2,1)
[~,sid]=sort(DEPTH);
imagesc(ACG_ALL(sid,:),[0 0.4])
title('V1')
colormap summer

ACG_ALL=[];
DEPTH = [];

for iF=1:size(rs_files,1)
    data = load(fullfile(savepath,rs_files{iF,1}{1}));
    data = data.data_out;
    raw_data = load(fullfile(raw_data_dir,rs_files{iF}{1}));
    
    wft =raw_data.sp.waveform_metrics;
    valid_rows = ismember(wft.cluster_id,data.good_cells);
    chan_number = wft.peak_channel(valid_rows);
    
    idx = ismember(chan_number,rs_files{iF,2}{1});
    %[~,sid]=sort(data.chan_number(idx),'descend');
     tmp = chan_number(idx);
    tmp = tmp-min(tmp);
    DEPTH = cat(1,DEPTH,tmp);
    ACG = data.xcorrs(idx,:);
    %ACG = ACG(sid,:);
    ACG_ALL=cat(1,ACG_ALL,ACG);
end
subplot(1,2,2)
[~,sid]=sort(DEPTH);
imagesc(ACG_ALL(sid,:),[0 0.4])
ma = max(size(ACG_ALL,1),ma);

title('RSC')
colormap summer
for ii=1:2
    set(subplot(1,2,ii),'YLim',[0 ma])
end
%saveas(gcf,'C:\Users\giocomolab\Desktop\heatmaps.pdf')