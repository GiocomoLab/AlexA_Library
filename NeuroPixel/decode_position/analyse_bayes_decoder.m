matfiles = dir('F:\attialex\tbtxcorr_decoder_0.5_Bayes_05/*.mat');

errorMat = [];
distanceMat = [];
error_time_mat = [];

take_idx_time = -450:450;

x_vec = 1:2:399;
ops = load_default_opt;
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.xbinedges = ops.edges;
ops.xbincent = .5*ops.edges(1:end-1)+.5*ops.edges(2:end);
savefig  = false;

REGIONS = {};
PARAMS = [];
for iF=1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name),'speed','trials','time_errorUniform','region');
    
    %errorMat = cat(3,errorMat,squeeze(data_out.scoreMat(1,:,:)));
   
    t1=data_out.time_errorUniform*-1;
    trials = data_out.trials-(data_out.trials(1))+1;
    bl_idx = trials<7 & abs(t1)'<350 & data_out.speed>5;
    a=polyfit(data_out.speed(bl_idx),t1(bl_idx)',1);
    REGIONS = cat(1,REGIONS,data_out.region(1));
    PARAMS = cat(1,PARAMS,a);
    
    gain_trial_onset = strfind(trials'==7,[0 1]);
    
    
    error_time_mat = cat(1,error_time_mat,t1(gain_trial_onset+take_idx_time));
    
    
    
end
%%
figure
histogram(PARAMS(startsWith(REGIONS,'MEC'),1))
hold on
histogram(PARAMS(startsWith(REGIONS,'VISp'),1))
