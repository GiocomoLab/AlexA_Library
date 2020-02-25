

ops.factors = -.25:0.01:.25;
ops.edges = 0:2:400;
ops.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = 3:20;
ops.TimeBin = 0.02;
ops.idx = [10:2:390]/2;% in bins
fi = gausswin(5);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = true;
OAK='/oak/stanford/groups/giocomo/';
%% savedir =
%savedir = fullfile(OAK,'attialex','speed_filtered_showMap');
savedir = fullfile('F:/temp/','speed_filtered');
imdir = fullfile(savedir,'images');
if ~isfolder(savedir)
    mkdir(savedir);
    
end
if ~isfolder(imdir)
    mkdir(imdir)
end
mf = matfile(fullfile(savedir,'parameters'),'Writable',true);
mf.ops = ops;


%% find files

% data_dir=fullfile(OAK,'attialex','NP_DATA');
% %session_name = {'AA5_190809_gain_1'};
% filenames = {};
% sn = dir(fullfile(data_dir,'*.mat'));
% for iS = 1:numel(sn)
%     if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
%         filenames{end+1}=sn(iS).name(1:end-4);
%     end
% end

gain = 0.8;
contrast = 100;
region = 'VISp';
[filenames,triggers] = getFilesCriteria(region,contrast,gain,'F:/NP_DATA');

%%
p=gcp('nocreate');
if isempty(p)
    parpool(6);
end

%%


for iF=1:numel(filenames)
    try
        data = load(filenames{iF});
        ops_temp = ops;
        ops_temp.trials=triggers{iF}(iRep)+[-18:3];
        test_trials = triggers{iF}(iRep)+[4:13];
        %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
        [data_out,fighandles] = findBestShifts(data,ops_temp);
        [~,mi]=max(data_out.all_stability,[],2);
        factors = ops_temp.factors(mi);
        
        
        
        [~,session_name,~]=fileparts(filenames{iF});
        if ops.plotfig
            for ifig = 1:numel(fighandles)
                saveas(fighandles{ifig},fullfile(imdir,sprintf('%s_%d.png',session_name,ifig)))
                close(fighandles{ifig})
            end
        end
        mf =matfile(fullfile(savedir,session_name));
        mf.stability = data_out.all_stability;
        mf.region = data_out.region;
        mf.test_trials = test_trials;
        mf.correlation_noshift = correlation_noshift;
        mf.correlation_shifted = correlation_shifted;
    catch ME
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
    end
end
