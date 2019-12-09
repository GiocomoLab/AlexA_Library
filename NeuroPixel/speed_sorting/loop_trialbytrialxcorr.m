data_dir=fullfile('F:','NP_DATA');
%session_name = {'AA5_190809_gain_1'};
session_name = {};
sn = dir(fullfile(data_dir,'*.mat'));
for iS = 1:numel(sn)
    if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
        session_name{end+1}=sn(iS).name(1:end-4);
    end
end

save_dir = fullfile('F:','temp','xcorrSpeed2');
if ~isfolder(save_dir)
    mkdir(save_dir)
end
%%

for iF=1:numel(session_name)
    data = load(fullfile(data_dir,session_name{iF}));
    try
    data_out = calculateTrialByTrialXCorr(data,5:20);
    if ~isempty(fieldnames(data_out))
        save(fullfile(save_dir,session_name{iF}),'data_out')
    end
    catch
    end
end