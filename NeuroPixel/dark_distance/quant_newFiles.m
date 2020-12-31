data_dir = 'F:\Alex\distance_tuning_xcorr_only';
raw_data_dir = 'F:\Alex\matfiles_new';
files = dir(fullfile(data_dir,'*.mat'));


%%
visp_files={    
   % {'AA_200919_1_MMdark_201010_13-32-24.mat'},{80:160};
    {'AA_200919_1_MMdark_201011_16-11-31.mat'},{85:200};
    {'AA_200919_3_MMdark_201019_12-00-55.mat'},{70:150};
    {'AA_200920_5_MMdark_201015_16-59-26.mat'},{85:180};
    };
 %%
 rs_files = {
    {'AA_200919_1_MMdark_201013_12-03-20.mat'},{130:250};
    %{'AA_200919_1_MMdark_201014_12-29-09.mat'},{};
    {'AA_200919_3_MMdark_201015_13-47-57.mat'},{150:300};
    {'AA_200920_5_MMdark_201018_11-48-04.mat'},{150:300};
    %{'AA_200920_5_MMdark_201019_14-21-28.mat'},
    };
%%
files = visp_files;
hmfig = figure()
tracefig = figure()
ops = load_mismatch_opt;
ops.dark = true;

for iF=1:size(files,1)
    data = load(fullfile(data_dir,files{iF}{1}));
    [~,sn]=fileparts(files{iF}{1});
    raw_data = load(fullfile(raw_data_dir,files{iF}{1}));
    good_cells = raw_data.sp.cids(raw_data.sp.cgs==2);
    frMat = calcFRVsDist(good_cells,[],raw_data,ops);
    
    [~,sid]=sort(data.chan_number,'descend');
    frMat=frMat(sid,:);
    ACG = data.ACG(sid,:);
    
    
    PI=[];
    figure(hmfig)
    subplot(1,2,1)
    imagesc(ACG,[0 0.4])
    hold on
    tmp = [];
    tmp_ub=[];
    tmp_fr = [];
    upper_limits = squeeze(data.acg_upper_limits(sid,:,5));
    firing_rate = data.firing_rate(sid);
    for ii=1:size(ACG,1)
        upper_bound = upper_limits(ii,:);
        [p,ip]=findpeaks(ACG(ii,:),'MinPeakProminence',.1,'SortStr','descend');
        if ~isempty(p)
            PI(end+1)=ip(1);
            
            if  p(1)>upper_bound(ip(1)) && firing_rate(ii)>1
                plot(ip(1),ii,'kx')
                tmp = cat(1,tmp,ACG(ii,:));
                tmp_ub = cat(1,tmp_ub,upper_bound);
                tmp_fr = cat(1,tmp_fr,frMat(ii,:));
            else
                plot(ip(1),ii,'rx')
            end
        end
    end
    
    colormap summer
    
    subplot(1,2,2)
    imagesc(tmp,[0 0.4])
    ylabel('ventral -> dorsal')
    title(sn,'Interpreter','None')
    figure(tracefig)
    for iT = 1:min(size(tmp,1),4)
        subplot(8,1,iT)
        plot(0:2:800,tmp(iT,:))
        hold on
        plot(0:2:800,tmp_ub(iT,:))
        subplot(8,1,iT+4)
        plot(tmp_fr(iT,:))
    end
    subplot(8,1,iT)
    legend({'Autocorrelation','99th prct shuffle'})
%scatter(raw_data.sp.waveform_metrics.amplitude,raw_data.sp.waveform_metrics.peak_channel)

    pause
    figure(tracefig)
    clf
    figure(hmfig)
    clf
end
%%

files = visp_files;
hmfig = figure()
ops = load_mismatch_opt;
ops.dark = true;

for iF=1:size(files,1)
    data = load(fullfile(data_dir,files{iF}{1}));
    [~,sn]=fileparts(files{iF}{1});
    
    IDX =ismember(data.chan_number,files{iF,2}{1});
    [~,sid]=sort(data.chan_number(IDX),'descend');
    ACG = data.ACG(IDX,:);
    ACG = ACG(sid,:);
    
    PI=[];
    subplot(1,3,iF)
    tmp = [];
    tmp_ub=[];
    tmp_fr = [];
    upper_limits =squeeze(data.acg_upper_limits(IDX,:,5));
    upper_limits = upper_limits(sid,:);
    for ii=1:size(ACG,1)
        upper_bound = upper_limits(ii,:);
        [p,ip]=findpeaks(ACG(ii,:),'MinPeakProminence',.2,'SortStr','descend');
        if ~isempty(p)
            PI(end+1)=ip(1);
            
            if  p(1)>upper_bound(ip(1))
               
                tmp = cat(1,tmp,ACG(ii,:));
                tmp_ub = cat(1,tmp_ub,upper_bound);
            
            end
        end
    end
    imagesc(tmp,[0 0.4])

    ylim([0,11])
    title(sprintf('%d out of %d',size(tmp,1),nnz(IDX)))
    colormap summer
end

