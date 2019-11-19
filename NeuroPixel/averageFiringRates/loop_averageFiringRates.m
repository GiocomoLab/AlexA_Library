%get files
region= 'MEC';
contrast = 100;
gain_to_look_at=1;
%[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/users/attialex/Desktop/data');
[filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,'/oak/stanford/groups/giocomo/attialex/NP_DATA');


savepath_root = '/oak/stanford/groups/giocomo/attialex/Images/avFiringRate';
%savepath_root = '/users/attialex/tmp/';
savepath = fullfile(savepath_root,sprintf('%s_%.2f_%d',region,gain_to_look_at,contrast));
if ~isfolder(savepath)
    mkdir(savepath)
end
p=gcp('nocreate');
if isempty(p)
    parpool(12);
end
%%
FiringRates = nan(numel(filenames),200);
FiringRatesNormalized = FiringRates;
FiringRatesIncreasing=FiringRates;
Names = cell(numel(filenames),1);
NUnits = nan(numel(filenames),1);
parfor iF=1:numel(filenames)
    data=load(filenames{iF});
    [tmp,tmpN,nunits,respIncreasing] = getSpatialMap(data,region);
    FiringRates(iF,:)=tmp;
    FiringRatesNormalized(iF,:)=tmpN;
    FiringRatesIncreasing(iF,:)=respIncreasing;
    [~,tmp,~]=fileparts(filenames{iF});
    Names{iF}=tmp;
    NUnits(iF)=nunits;
end
    
out.FiringRates = FiringRates;
out.FiringRatesNorm = FiringRatesNormalized;
out.FiringRatesIncreasing = FiringRatesIncreasing;
out.names = Names;
out.NUnits= NUnits;
save(fullfile(savepath,'firingRateData.mat'),'out')
%%
figure
x=1:2:400;

idx = (1:199);
subplot(2,1,1)
errorbar(x(idx),mean(out.FiringRates(:,idx)*50),std(out.FiringRates(:,idx)*50)/sqrt(numel(out.NUnits)))
hold on
errorbar(x(idx),nanmean(out.FiringRatesIncreasing(:,idx)*50),nanstd(out.FiringRatesIncreasing(:,idx)*50)/sqrt(numel(out.NUnits)));
xlabel('Position [cm]')
ylabel('Firing Rate [Hz]')
subplot(2,1,2)
errorbar(x(idx),mean(out.FiringRatesNorm(:,idx)),std(out.FiringRatesNorm(:,idx))/sqrt(numel(out.NUnits)))
xlabel('Position [cm]')
ylabel('Z Scored Rate')
