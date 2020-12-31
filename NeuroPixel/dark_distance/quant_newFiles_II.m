data_dir = 'F:\Alex\distance_tuning_xcorr_only_3';
raw_data_dir = 'F:\Alex\matfiles_new';
files = dir(fullfile(data_dir,'*.mat'));


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
files = visp_files;
hmfig = figure();
tracefig = figure();
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
    imagesc([0:2:800],1:size(ACG,1),ACG,[0 0.4])
    hold on
    tmp = [];
    tmp_ub=[];
    tmp_fr = [];
    upper_limits = squeeze(data.acg_upper_limits(sid,:,5));
    firing_rate = data.firing_rate(sid);
    for ii=1:size(ACG,1)
        upper_bound = upper_limits(ii,:);
        %[p,ip]=findpeaks(ACG(ii,:),'MinPeakProminence',.1,'SortStr','descend');
        if ~isnan(data.peak_all(sid(ii)))
       
            if  data.peak_prom_all(sid(ii))>0.1 && data.pval(sid(ii))<=0.01
                plot(data.peak_loc_all(sid(ii)),sid(ii),'kx')
                tmp = cat(1,tmp,ACG(ii,:));
                tmp_ub = cat(1,tmp_ub,upper_bound);
                tmp_fr = cat(1,tmp_fr,frMat(ii,:));
            else
                plot(data.peak_loc_all(sid(ii)),sid(ii),'rx')
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
    if ~isempty(iT)
    subplot(8,1,iT)
        legend({'Autocorrelation','99th prct shuffle'})

    end
%scatter(raw_data.sp.waveform_metrics.amplitude,raw_data.sp.waveform_metrics.peak_channel)

    pause
    figure(tracefig)
    clf
    figure(hmfig)
    clf
end
%%

files = visp_files;
hmfig = figure();
tracefig = figure();
ops = load_mismatch_opt;
ops.dark = true;
PROM_ALL=[];
PEAK_ALL=[];
TUNED_ALL=[];
PEAK_SHUFF_ALL=[];
for iF=1:size(files,1)
    data = load(fullfile(data_dir,files{iF}{1}));
    
    idx = ismember(data.chan_number,files{iF,2}{1});
    PROM_ALL=cat(1,PROM_ALL,data.peak_prom_all(idx));
    PEAK_ALL = cat(1,PEAK_ALL,data.peak_all(idx));
    tmp = data.peak_prom_all>0.1 & data.pval<0.01;
    TUNED_ALL=cat(1,TUNED_ALL,tmp(idx));
    PEAK_SHUFF_ALL=cat(1,PEAK_SHUFF_ALL,data.peak_shuf(idx,:));
    
end
peak_zscore=nan(size(PEAK_ALL));
for iC=1:size(PEAK_ALL,1)
    peak_zscore(iC)=(PEAK_ALL(iC)-nanmean(PEAK_SHUFF_ALL(iC,:)))/nanstd(PEAK_SHUFF_ALL(iC,:));
end

figure
hold on
scatter(PROM_ALL(TUNED_ALL==1),peak_zscore(TUNED_ALL==1),45,'b','.')
scatter(PROM_ALL(TUNED_ALL==0),peak_zscore(TUNED_ALL==0),45,'k','.')
fprintf('%d out of %d tuned \n',nnz(TUNED_ALL),numel(TUNED_ALL));


%%
figure
ACG_ALL=[];
DEPTH = [];
ma = 0;
for iF=1:size(visp_files,1)
    data = load(fullfile(data_dir,visp_files{iF,1}{1}));
    idx = ismember(data.chan_number,visp_files{iF,2}{1});
    %[~,sid]=sort(data.chan_number(idx),'descend');
    tmp = data.chan_number(idx);
    tmp = tmp-min(tmp);
    DEPTH = cat(1,DEPTH,tmp);
    ACG = data.ACG(idx,:);
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
    data = load(fullfile(data_dir,rs_files{iF,1}{1}));
    idx = ismember(data.chan_number,rs_files{iF,2}{1});
    %[~,sid]=sort(data.chan_number(idx),'descend');
     tmp = data.chan_number(idx);
    tmp = tmp-min(tmp);
    DEPTH = cat(1,DEPTH,tmp);
    ACG = data.ACG(idx,:);
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
saveas(gcf,'C:\Users\giocomolab\Desktop\heatmaps.pdf')