%waveform extractor

%% get files


%%

for iF =1:numel(filenames)
    data = load(filenames{iF});
    dat_path = data.sp.dat_path;
    %find full path
    fullpath = '';
    
    [mean_waveforms,mean_temlate_waveforms]=get_waveforms(realpath);
    
    %data.mean_waveforms = mean_waveforms;
    %data.mean_template_waveforms = mean_template_waveforms;
    save(filenames{iF},'mean_waveforms,'-append')
    
end