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
Names = cell(numel(filenames),1);
NUnits = nan(numel(filenames),1);
parfor iF=1:numel(filenames)
    data=load(filenames{iF});
    [tmp,tmpN,nunits] = getSpatialMap(data,region);
    FiringRates(iF,:)=tmp;
    FiringRatesNormalized(iF,:)=tmpN;
    [~,tmp,~]=fileparts(filenames{iF});
    Names{iF}=tmp;
    NUnits(iF)=nunits;
end
    
out.FiringRates = FiringRates;
out.FiringRatesNorm = FiringRatesNormalized;
out.names = Names;
out.NUnits= NUnits;
save(fullfile(savepath,'firingRateData.mat'),'out')