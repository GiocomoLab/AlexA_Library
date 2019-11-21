
%get files
region= 'MEC';
contrast = 100;
gain_to_look_at=0.5;
%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/users/attialex/Desktop/data');
[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
%%
savepath = '/oak/stanford/groups/giocomo/attialex/Images/xcorrPrePostGain';
if ~isfolder(savepath)
    mkdir(savepath)
end
%load data
for iF=1:numel(filenames)
    data = load(filenames{iF});
    trials = triggers{iF}(1):triggers{iF}(1)+3;
    trials_xtx = (triggers{iF}(1)-4):triggers{iF}(1)+3;
    [ACG_corr_gain,ZERO_LAG_PRE,ZERO_LAG_POST] = extractXCorr(data,trials,region);
    [ACG_corr_BL,ZERO_LAG_PRE_BL,ZERO_LAG_POST_BL] = extractXCorr(data,trials-4,region);
    XTX_GAIN = getCovMatrix(data,region,trials_xtx,2,1:4,-1);
    XTX_BL = getCovMatrix(data,region,trials_xtx-4,2);
    h=figure();
    subplot(2,3,1)
    scatter(ACG_corr_gain,ACG_corr_BL)
    axis image
    xlim([-.2 1])
    ylim([-.2 1])
    xlabel('ACG_gain')
    ylabel('ACG_BL')
    title('Similarity of ACG')
    set(gca,'YTick',[0 0.5 1])
    set(gca,'YTick',[0 0.5 1])

    grid on
    subplot(2,3,2)
    scatter(ZERO_LAG_PRE,ZERO_LAG_POST)
    xlabel('Zero lag Pre')
    ylabel('Zero lag Post')
    title('Zero lag Gain')
        axis image
    ma=max([ZERO_LAG_PRE;ZERO_LAG_POST]);
    axis([0 ma 0 ma])
    grid on
    subplot(2,3,3)
    scatter(ZERO_LAG_PRE_BL,ZERO_LAG_POST_BL)
    axis image
    ma=max([ZERO_LAG_PRE_BL;ZERO_LAG_POST_BL]);
    axis([0 ma 0 ma])
        xlabel('Zero lag Pre')
    ylabel('Zero lag Post')
    grid on
    title('Zero Lag Baseline')
    
    subplot(2,3,5)
    imagesc(XTX_GAIN,[0 1])
    axis image
    subplot(2,3,6)
    imagesc(XTX_BL,[0 1])
    axis image
    [~,session_name,~]=fileparts(filenames{iF});
    
    saveas(h,fullfile(savepath,sprintf('%s.png',session_name)))
    close(h)
end
%extract good cells


