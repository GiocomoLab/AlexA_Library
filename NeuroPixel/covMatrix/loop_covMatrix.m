files = dir('F:\NP_DATA\np*_gain*.mat');

nF = numel(files);
trials = [1:30];
nTrials = numel(trials);
XTX = nan(nBins*nTrials,nBins*nTrials,nF);
GAINS = zeros(nTrials,nF);
CONTRASTS = zeros(nTrials,nF);
%parpool(4);

for iF=1:length(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    
    session_name = files(iF).name(1:end-4);
    this_sess = strcmp(cell_info.Session,session_name);
    
    [tmpXTX,trial,gain,contrast]=getCovMatrix(data,cell_info(this_sess,:),trials);
    XTX(:,:,iF)=tmpXTX;
    GAINS(:,iF)=gain;
    CONTRASTS(:,iF)=contrast;
end
    