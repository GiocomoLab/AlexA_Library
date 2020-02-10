
ops.factors = -.25:0.01:.25;
ops.edges = 0:2:400;
ops.search_range=[11:189]; % in bins

ops.midpoints = ops.edges(1:end-1)+(ops.edges(1)+ops.edges(2))/2;
ops.BinWidth = 2;
ops.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = 3:24;
ops.TimeBin = 0.02;
ops.speedWindow = -10:-1; % in cm

fi = gausswin(3);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = true;
OAK='/oak/stanford/groups/giocomo/';

%%
% data_dir=fullfile('F:','NP_DATA');
% %session_name = {'AA5_190809_gain_1'};
% filenames = {};
% sn = dir(fullfile(data_dir,'*.mat'));
% for iS = 1:numel(sn)
%     if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
%         filenames{end+1}=sn(iS).name(1:end-4);
%     end
% end
[filenames,triggers] = getFilesCriteria('VISp',100,0.8,'F:/NP_DATA');


%%
output = cell(1,numel(filenames));
for iF=1:numel(filenames)
    try
    %data = load(fullfile(data_dir,filenames{iF}));
    data = load(filenames{iF});
    ops.trials=triggers{iF}(1)+[-18:3];
    data_out = findPeakAlignement(data,ops);
    
    output{iF}=data_out;
    catch
        disp('whoopsie')
    end
end

%%

FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
reg= [];
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    if numel(output{iF}.stability) ~= numel(output{iF}.region)
        disp(iF)
        continue
    end
    FAST =cat(1,FAST,output{iF}.all_fast);
    SLOW = cat(1,SLOW,output{iF}.all_slow);
    GAIN = cat(1,GAIN,output{iF}.all_gain);
    STAB = cat(1,STAB,output{iF}.stability);
    reg = cat(2,reg,output{iF}.region);
end


figure
idx = STAB>.2 & ismember(reg,'VISp')';
plot(nanmean(FAST(idx,:)))
hold on
plot(nanmean(SLOW(idx,:)));
plot(nanmean(GAIN(idx,:)))
legend({'fast','slow','gain'})
%%
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
reg= [];
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    if numel(output{iF}.stability) ~= numel(output{iF}.region)
        disp(iF)
        continue
    end
    FAST =cat(1,FAST,output{iF}.all_fast);
    SLOW = cat(1,SLOW,output{iF}.all_slow);
    GAIN = cat(1,GAIN,output{iF}.all_gain);
    STAB = cat(1,STAB,output{iF}.stability);
    reg = cat(2,reg,output{iF}.subregion);
end


figure
idx = STAB>.1 & ismember(reg,'VISp2/3')';
plot(nanmean(FAST(idx,:)))
hold on
plot(nanmean(SLOW(idx,:)));
plot(nanmean(GAIN(idx,:)))
legend({'fast','slow','gain'})