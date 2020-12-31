%%
cl_filenames=dir('Z:\giocomo\attialex\NP_DATA_2\*_mismatch*.mat');

%%
% n=0;
% fi = {};
% for iF=1:numel(cl_filenames)
%     fn_new = strrep(cl_filenames(iF).name,'mismatch','playback');
%     if isfile(fullfile(cl_filenames(1).folder,fn_new))
%         n=n+1;
%         fi{end+1}=fn_new;
%     end
% end

%%
opt = load_mismatch_opt;

for iF=1:size(cl_filenames,1)
    fn_new = strrep(cl_filenames(iF).name,'mismatch','playback');
    if ~isfile(fullfile(cl_filenames(iF).folder,fn_new))
        continue
    end
    data_cl = load(fullfile(cl_filenames(iF).folder,cl_filenames(iF).name));
    data_ol = load(fullfile(cl_filenames(iF).folder,fn_new));
    good_cells = data_ol.sp.cids(data_ol.sp.cgs==2);
    try
        glmDataPB = fitGLM_OLCL_oldDataFormat(data_ol,data_cl,good_cells,1);
        
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
        
        p=cat(2,glmDataPB(:).yhat_cl);
        snps = extract_snps(p',mm_trigs,'win',[-100 149]);
        
        run_ons = strfind(smooth_speed>opt.speed_t,[zeros(1,30),ones(1,30)])+30;
        [spikeTimes,~,aux,~,count_vec_run]=extract_triggered_spikeTimes(data.sp,data.post(run_ons),'cluIDs',good_cells,'win',opt.extract_win,'aux',[data.post' ;smooth_speed],'aux_win',opt.aux_win);
        avg_yhat = squeeze(mean(snps,3));
        if isfield(data,'anatomy')
            if isfield(data.anatomy,'parent_shifted')
                reg = data.anatomy.parent_shifted;
            else
                reg = data.anatomy.cluster_parent;
            end
            if iscolumn(reg)
                reg = reg';
            end
            
            reg=reg(data.sp.cgs==2);
        else
            reg = {};
        end
        mf = matfile(['F:\Alex\glmFits_oldData2\' fn_new],'Writable',true);
        mf.glmData = glmDataPB;
        mf.mm_resp = count_vec;
        mf.mm_predicted = avg_yhat;
        mf.reg = reg;
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

