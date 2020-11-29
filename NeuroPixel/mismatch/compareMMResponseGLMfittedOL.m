matfiles = dir('F:\Alex\matfiles_new\*_PB_*.mat');

for iF = 1:numel(matfiles)
    pb_file = matfiles(iF).name;
    dashes = strfind(pb_file,'_');
    pb_files{iF,1}=pb_file;
    
    if numel(dashes)<5
        continue
    end
    
    animal = pb_file(1:dashes(3)-1);
    
    mm_string = pb_file(1:dashes(5));
    mm_string = strrep(mm_string,'_PB_','_MM_');
    pot = dir(fullfile(matfiles(iF).folder,[mm_string '*']))
    
    pb_files{iF,2} = pot(1).name;
end

%%
data_root = 'F:\Alex\matfiles_new'
files={
    {'AA_200919_1_PB_11-01-40.mat'       },    {'AA_200919_1_MM_10-47-17.mat'                          };
    {'AA_200919_1_PB_18-03-39.mat'       },    {'AA_200919_1_MM_17-41-26.mat'                          };
    {'AA_200919_1_PB_201007_15-42-13.mat'},    {'AA_200919_1_MM_201007_15-28-09.mat'};
    {'AA_200919_1_PB_201010_13-06-38.mat'},    {'AA_200919_1_MM_201010_12-39-29.mat'};
    {'AA_200919_1_PB_201011_16-45-13.mat'},    {'AA_200919_1_MM_201011_16-24-35.mat'};
    {'AA_200919_1_PB_201013_12-41-25.mat'},    {'AA_200919_1_MM_201013_12-16-45.mat'};
    {'AA_200919_1_PB_201014_12-57-18.mat'},    {'AA_200919_1_MM_201014_12-41-44.mat'};
    {'AA_200919_3_PB_201010_14-58-16.mat'},    {'AA_200919_3_MM_201010_14-43-47.mat'};
    {'AA_200919_3_PB_201011_14-40-18.mat'},    {'AA_200919_3_MM_201011_14-25-36.mat'};
    {'AA_200919_3_PB_201013_16-49-24.mat'},    {'AA_200919_3_MM_201013_16-36-05.mat'};
    {'AA_200919_3_PB_201015_13-33-47.mat'},    {'AA_200919_3_MM_201015_13-19-37.mat'};
    {'AA_200919_3_PB_201018_13-22-13.mat'},    {'AA_200919_3_MM_201018_13-06-15.mat'};
    {'AA_200919_3_PB_201019_11-47-20.mat'},    {'AA_200919_3_MM_201019_11-33-43.mat'};
    {'AA_200920_4_PB_13-47-38.mat'       },    {'AA_200920_4_MM_13-14-15.mat'                          };
    {'AA_200920_5_PB_201010_17-05-15.mat'},    {'AA_200920_5_MM_201010_16-46-40.mat'};
    {'AA_200920_5_PB_201011_12-45-43.mat'},    {'AA_200920_5_MM_201011_12-22-36.mat'};
    {'AA_200920_5_PB_201013_14-55-16.mat'},    {'AA_200920_5_MM_201013_14-33-40.mat'};
    {'AA_200920_5_PB_201014_14-34-42.mat'},    {'AA_200920_5_MM_201014_14-21-14.mat'};
    {'AA_200920_5_PB_201015_16-44-28.mat'},    {'AA_200920_5_MM_201015_16-29-02.mat'};
    {'AA_200920_5_PB_201018_11-34-12.mat'},    {'AA_200920_5_MM_201018_11-20-13.mat'};
    {'AA_200920_5_PB_201019_14-03-48.mat'},    {'AA_200920_5_MM_201019_13-46-40.mat'};
    };


%%
    opt = load_mismatch_opt;

for iF=1:size(files,1)
    
    data_cl = load(fullfile(data_root,files{iF,2}{1}));
    data_ol = load(fullfile(data_root,files{iF,1}{1}));
    good_cells = data_ol.sp.cids(data_ol.sp.cgs==2);
    try
    glmDataPB = fitGLM_OLCL(data_ol,data_cl,good_cells,1);
    
    data = data_cl;
    if isfield(data,'vr_data_resampled')
        mismatch_trigger = data.vr_data_resampled.MM>0.5;
        true_speed = data.vr_data_resampled.velM;
    else
        mismatch_trigger = data.mismatch_trigger;
        mismatch_trigger = mismatch_trigger==0;
        true_speed = data.true_speed;
    end
    if iscolumn(mismatch_trigger)
        mismatch_trigger = mismatch_trigger';
    end
    %good_cells = data.sp.ks_cluster.cluster_id(startsWith(data.sp.rf_cluster.group,'good'));
    %%
    opt.speed_t=0.05;
    all_mm_trigs=strfind(mismatch_trigger>0.9,[0 0 1 1])+2;
    if iscolumn(true_speed)
        speed=true_speed';
    else
        speed=true_speed;
    end
    filt = gausswin(61); %61 pretty close to what we use in other
    filt = filt/sum(filt);
    smooth_speed = conv(speed,filt,'same');
    %%
    run_periods=smooth_speed>opt.speed_t;
    run_window=-30:30;
    possibles=strfind(run_periods,ones(1,length(run_window)))+floor(.5*length(run_window));
    aux_vec = [data.post' ;smooth_speed];

    mm_trigs=all_mm_trigs(ismember(all_mm_trigs,possibles));
    %%

    [spikeTimes,~,aux,~,count_vec]=extract_triggered_spikeTimes(data.sp,data.post(mm_trigs),'cluIDs',good_cells,'win',opt.extract_win,'aux',aux_vec,'aux_win',opt.aux_win);
    
    p=cat(2,glmDataPB(:).yhat);
    snps = extract_snps(p',mm_trigs,'win',[-100 149]);
    
    run_ons = strfind(smooth_speed>opt.speed_t,[zeros(1,30),ones(1,30)])+30;
    [spikeTimes,~,aux,~,count_vec_run]=extract_triggered_spikeTimes(data.sp,data.post(run_ons),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed],'aux_win',opt.aux_win);
    avg_yhat = squeeze(mean(snps,3));
   
    mf = matfile(['F:\Alex\glmFits\' files{iF,2}{1}],'Writable',true);
    mf.glmData = glmDataPB;
    mf.mm_resp = count_vec;
    mf.mm_predicted = avg_yhat;
    figure;
    plot(opt.time_vecs,mean(count_vec))
    hold on
    plot(opt.time_vecs,mean(avg_yhat)/0.02)
    title(iF)
    catch ME
        disp(ME.message)
    end
end

%%

