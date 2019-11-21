%get files
region= 'MEC';
contrast = 100;
gain_to_look_at=1;
%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/users/attialex/Desktop/data');
[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/oak/stanford/groups/giocomo/attialex/NP_DATA');


%parpool(4);
p=gcp('nocreate');
if isempty(p)
    parpool(12);
end
%%
nF = numel(filenames);
trials = [3:20];
nTrials = numel(trials);
nBins = 200;
XTX = nan(nBins*nTrials,nBins*nTrials,nF);
GAINS = zeros(nTrials,nF);
CONTRASTS = zeros(nTrials,nF);
NAMES  = cell(numel(filenames),1);


parfor iF=1:length(filenames)
    data = load(filenames{iF});
    
    [~,session_name,~] = fileparts(filenames{iF});
 
    
    [tmpXTX,trial,gain,contrast,nunits]=getCovMatrix(data,region,trials,2);
    XTX(:,:,iF)=tmpXTX;
    GAINS(:,iF)=gain;
    CONTRASTS(:,iF)=contrast;
    NAMES{iF}=session_name;
    NUnits(iF)=nunits;
end
%%

tmp = nanmean(XTX,3);
plot(diag(tmp,200))
%%
AVG=zeros(200,size(XTX,3));
for iS =1:size(XTX,3)
    tmp = diag(XTX(:,:,iS),-200);
    tmpM = reshape(tmp,200,17);
    AVG(:,iS)=mean(tmpM,2);
   
end
%%
set(0,'DefaultFigureRenderer','Painters')
tmpA=mean(AVG,2);
tmpA=circshift(tmpA,100);
tmpS=std(AVG,[],2);
tmpS=tmpS/sqrt(size(XTX,3));
tmpS=circshift(tmpS,100);
figure
boundedline(-198:2:200,tmpA,tmpS,'alpha')
ylabel('Similarity')
xlabel('Distance from Reward Location [cm]')
saveas(gcf,'/oak/stanford/groups/giocomo/attialex/Images/xcov_shifted.pdf')

%%
set(0,'DefaultFigureRenderer','Painters')
tmpA=mean(AVG,2);
tmpS=std(AVG,[],2);
tmpS=tmpS./sqrt(size(XTX,3));
figure
boundedline(2:2:400,tmpA,tmpS,'alpha')
ylabel('Similarity')
xlabel('Track Position [cm]')
%%
names = {};
for iF=1:numel(filenames)
    [~,name,~]=fileparts(filenames{iF});
    names{iF}=name(1:4);
end
nnames = numel(unique(names));
avunits = mean(NUnits);
sprintf('%d Mice, %d Sessions, %.2f cells',nnames,numel(filenames),avunits)