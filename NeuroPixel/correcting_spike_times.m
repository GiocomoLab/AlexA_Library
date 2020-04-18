NP_DIR = 'F:/NP_DATA';
matfiles = dir(fullfile(NP_DIR,'npJ3*.mat'))
dest_path = 'F:/NP_DATA2';
if ~isfolder(dest_path)
    mkdir(dest_path)
end
malcolm_dir = 'Z:\giocomo\export\data\Projects\ContrastExperiment_neuropixels'
test_name = 'npJ5_0504_baseline_gain_3.mat'
%%
for iF=1:numel(matfiles)
    sn = matfiles(iF).name;
    [~,session_name]=fileparts(matfiles(iF).name);
    parts = strsplit(session_name,'_');
    
    animal = parts{1};
    animal = animal(3:end);
    %find it on oak
    %position file:
    pos_file = '';
    for ip=2:numel(parts)
        pos_file = strcat(pos_file,parts{ip},'_');
    end
    pos_file = [pos_file 'position.txt'];
    
    
    fi = dir(fullfile(malcolm_dir,animal,'neuropixels_data','*',sn));
    
    if numel(fi)==1
        data_dir = fi.folder;
        [~,main_name]=fileparts(data_dir);
        
        NIDAQ_file = dir(fullfile(data_dir,'*nidq.bin'));
        NIDAQ_file = fullfile(data_dir,NIDAQ_file(1).name);
        NIDAQ_config = dir(fullfile(data_dir,'*nidq.meta'));
        NIDAQ_config = fullfile(data_dir,NIDAQ_config(1).name);
        spike_dir = fullfile(data_dir,strcat(main_name,'_imec0'));
        out = isfile(NIDAQ_file);
        sprintf('found %s \n, NIDAQ: %d',data_dir,out)
        
        
        %get the nidaq sample rate & get number of recorded nidaq channels
        dat=textscan(fopen(NIDAQ_config),'%s %s','Delimiter','=');
        names=dat{1};
        vals=dat{2};
        loc=contains(names,'niSampRate');
        sync_sampling_rate=str2double(vals{loc});
        
        loc2=contains(names,'nSavedChans');
        n_channels_nidaq=str2double(vals{loc2});
        
        % get neuropixels sync pulse times
        fpNIDAQ=fopen(NIDAQ_file);
        datNIDAQ=fread(fpNIDAQ,[n_channels_nidaq,Inf],'*int16');
        fclose(fpNIDAQ);
        syncDat=datNIDAQ(2,:)>1000;
        
        
        frame_times_np = find(abs(diff(syncDat))==1)+1;
        frame_times_np = frame_times_np/sync_sampling_rate;
        
        % read vr position data
        formatSpec = '%f%f%f%f%f%[^\n\r]';
        delimiter = '\t';
        if isfile(fullfile(data_dir,strcat(session_name,'_position.txt')))
            vr_session_name = session_name;
        elseif isfile(fullfile(data_dir,strcat(session_name(6:end),'_position.txt')))
            vr_session_name = session_name(6:end);
        end
        fn_vr = fullfile(data_dir,strcat(vr_session_name,'_position.txt'));
        fn_trial = fullfile(data_dir,strcat(vr_session_name,'_trial_times.txt'));
        fn_lick = fullfile(data_dir,strcat(vr_session_name,'_licks.txt'));
        
        %         else
        %             fn_vr = fullfile(data_dir,strcat(session_name,'_position.txt'));
        %             fn_trial = fullfile(data_dir,strcat(session_name,'_trial_times.txt'));
        %             fn_lick = fullfile(data_dir,strcat(session_name,'_licks.txt'));
        %
        %             %fn_vr = fullfile(data_dir,'..','..','VR',strcat(session_name,'_position.txt'));
        %             %fn_trial = fullfile(data_dir,'..','..','VR',strcat(session_name,'_trial_times.txt'));
        %             %fn_lick = fullfile(data_dir,'..','..','VR',strcat(session_name,'_licks.txt'));
        %         end
        fid = fopen(fn_vr,'r');
        dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
        fclose(fid);
        vr_position_data = cat(2,dataArray{1:5});
        %vr_position_data = vr_position_data(1:49334,:);
        nu_entries = nnz(~isnan(vr_position_data(1,:)));
        %vr_position_data=vr_position_data(49335:end,:);
        vr_ttl=vr_position_data(:,nu_entries); %assuming TTL in last and timestamp in second last column
        frame_times_vr=vr_position_data(:,nu_entries-1);
        
        % read vr trial data
        fid = fopen(fn_trial,'r');
        vr_trial_data = fscanf(fid, '%f', [4,inf])';
        fclose(fid);
        trial_contrast = [100; vr_trial_data(:,2)];
        trial_gain = [1; vr_trial_data(:,3)];
        num_trials = numel(trial_gain);
        
        % read vr licking data
        try
            fid = fopen(fn_lick,'r');
            vr_lick_data = fscanf(fid, '%f', [2,inf])';
            fclose(fid);
            lickx = vr_lick_data(:,1);
            lickt = vr_lick_data(:,2);
        catch ME
            fprintf('no lick file found \n')
            lickx=NaN;
            lickt=NaN;
        end
        
        
        % set vr frame times to be the time of neuropixels pulses
        % make sure the number of frames matches (can be off by one because of
        % odd/even numbers of frames)
        
        tmp_diff=diff(frame_times_np);
        [mm,step_idx]=find(tmp_diff>2); %%CHANGE BACK TO 2
        
        sess_length=diff([0 step_idx length(frame_times_np)]);
        midpoint = ([0 step_idx] + [step_idx length(frame_times_np)])/2;
        %step_idx=step_idx+1;
        frametimes_nlOld = frame_times_np;
        [~,ml]=min(abs(sess_length-numel(frame_times_vr)));
        if length(mm)>=1
            
            
            sess=ml;
            step_idx = [0 step_idx length(frame_times_np)];
            
            idx_start=step_idx(sess)+1;
            idx_stop = step_idx(sess+1);
            frame_times_np=frame_times_np(idx_start:idx_stop);
        else
            is_mismatch=0;
        end
        %%
        if abs(numel(frame_times_np) - numel(frame_times_vr)) <= 1
            idx=1:min(numel(frame_times_np),numel(frame_times_vr)); %use shorter index
            post = frame_times_np(idx)';
            vr_position_data=vr_position_data(idx,:);
            posx = vr_position_data(:,1);
            
            % transform lick times into neuropixels reference frame
            beta = [ones(size(post)) frame_times_vr(idx)]\post;
            lickt = beta(1) + lickt*beta(2);
        else
            disp('ERROR: number of sync pulses does not match number of frames.')
            continue
        end
        figure
        scatter(diff(post),diff(frame_times_vr(idx)),2,1:length(idx)-1)
        r=corrcoef(diff(post),diff(frame_times_vr(idx)));
        title(sprintf('corr coeff = %0.3f',r(1,2)));
        dest_pic = fullfile(dest_path,[session_name '.png']);
        saveas(gcf,dest_pic)
        close(gcf);
        offset = post(1);
        src_file = fullfile(matfiles(iF).folder,matfiles(iF).name);
        data_out = load(src_file);
        sp = data_out.sp;
        sp.vr_session_offset = offset;

        
        syncDatNIDAQ=datNIDAQ(1,:)>1000;
        % convert NIDAQ sync data into time data by dividing by the sampling rate
        ts_NIDAQ = strfind(syncDatNIDAQ,[0 1])/sync_sampling_rate;
        % ts_NIDAQ: these are the sync pulse times relative to the NIDAQ board
        % Now, we do the same, but from the perspective of the Imec board.
        LFP_config = dir(fullfile(spike_dir,'*.lf.meta'));
        LFP_config = fullfile(LFP_config.folder,LFP_config.name);
        LFP_file = dir(fullfile(spike_dir,'*.lf.bin'));
        LFP_file = fullfile(LFP_file.folder,LFP_file.name);
        dat=textscan(fopen(LFP_config),'%s %s','Delimiter','=');
        names=dat{1};
        vals=dat{2};
        loc=contains(names,'imSampRate');
        lfp_sampling_rate=str2double(vals{loc});
        % for loading only a portion of the LFP data
        fpLFP = fopen(LFP_file);
        fseek(fpLFP, 0, 'eof'); % go to end of file
        fpLFP_size = ftell(fpLFP); % report size of file
        fpLFP_size = fpLFP_size/(2*384);
        fclose(fpLFP);
        % get the sync pulse times relative to the Imec board
        fpLFP=fopen(LFP_file);
        fseek(fpLFP,384*2,0);
        ftell(fpLFP);
        datLFP=fread(fpLFP,[1,round(fpLFP_size/4)],'*int16',384*2); % this step used to take forever
        fclose(fpLFP);
        syncDatLFP=datLFP(1,:)>10;
        ts_LFP = strfind(syncDatLFP,[0 1])/lfp_sampling_rate;
        % ts_LFP: these are the sync pulse times relative to the Imec board
        % PART 2: TIME CORRECTION
        lfpNIDAQdif = ts_LFP - ts_NIDAQ(1:size(ts_LFP, 2)); % calculate the difference between the sync pulse times
        fit = polyfit(ts_LFP, lfpNIDAQdif, 1); % linear fit
        correction_slope = fit(1); % this is the amount of drift we get per pulse (that is, per second)
        % save the old, uncorrected data as sp.st_uncorrected and save the new,
        % corrected data as sp.st (as many of your analyses are using sp.st).
        
        %correct spike time
        
        sp.st_uncorrected = sp.st; % save uncorrected spike times (st)
        st_corrected = sp.st - (sp.st+offset) * correction_slope; % in two steps to avoid confusion
        sp.st = st_corrected; % overwrite the old sp.st
        valid_idx = sp.st>=0;
        sp.st=sp.st(valid_idx);
        sp.st_uncorrected=sp.st_uncorrected(valid_idx);
        sp.clu=sp.clu(valid_idx);
        sp.spikeTemplate = sp.spikeTemplates(valid_idx);
        sp.tempScalingAmps = sp.tempScalingAmps(valid_idx);
        dest_file = fullfile(dest_path,matfiles(iF).name);
        copyfile(src_file,dest_file);
        save(dest_file,'sp','-append');
    else
        %sprintf('did not find %s \n',pos_file)
    end
end
%'0504_baseline_gain_1_position'