%% struct_names
%first run loop_findshifts to generate the data
path = '/Volumes/Samsung_T5/attialex/speed_filtered_correctedData';
path = '/Volumes/T7/attialex/speed_filtered_correctedData_shortidx2';

filenames = dir(fullfile(path,'*.mat'));
ops = load(fullfile(path,'parameters.mat'));
ops = ops.ops;
opt = load_default_opt;
%% files with .8 gain (to make sure we only use files that were also used in gain=0.8 experiments
%% makes it easier afterwards
gain = 0.8;
contrast = 100;
regions = {'VISp','RS','MEC'};
filenames_08 = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile('/Volumes/T7/attialex/NP_DATA_corrected'));
    filenames_08=cat(2,filenames_08,tmp1);
    triggers = cat(2,triggers,tmp2);
end

NBlocks = [];
Frac_Good = [];
%%

STABILITY = struct();
GoodCells = struct();
animal_counter = 0;
ANIMALS = struct();
for iF=1:numel(filenames)
    try
        filename=filenames(iF).name;
        if nnz(endsWith(filenames_08,filename))==0
            continue
        end
        parts=strsplit(filename,'_');
        animal = parts{1};
        if ~isfield(ANIMALS,animal)
            animal_counter = animal_counter +1;
            ANIMALS.(animal) = animal_counter;
        end
        aid = ANIMALS.(animal);
        a=load(fullfile(path,filename));
        [unique_regions,~,ridx]=unique(a.region);
        NBlocks = cat(1,NBlocks,size(a.all_stability,1));
        stab = a.all_stability;
                stable_blocks = stab>=opt.stab_thresh;
                tmp_factors = a.all_factors;
                
                tmp_factors(~stable_blocks)=nan; %set shifts to nan where stability <.5
                tmp_factors2 = tmp_factors;
                valid_region = startsWith(a.region,{'MEC','RS','VISp'});
                p=a.all_stability(:,valid_region);
                q=nanmean(p,2);
                %f=sum(p>ops.stab_thresh,2);
                %tmp_factors2(f<55,:)=nan;
                tmp_factors2(q<0.5,:)=nan;
                enough_data_idx = sum(~isnan(tmp_factors),1)>2; %only use data where there are more than 2 stable blocks
                enough_data_idx2 = sum(~isnan(tmp_factors2),1)>2;
                tmp_idx = startsWith(a.region,regions);
        Frac_Good = cat(1,Frac_Good,[sum(~isnan(tmp_factors(:,tmp_idx)),1)', ones(nnz(tmp_idx),1)*size(a.all_stability,1)]);
        for iR=1:numel(unique_regions)
            
            r=unique_regions{iR};
            iidx = ridx==iR;
            if startsWith(r,'RSPagl')
                %keyboard
                r = 'RAGL';
            end
            if startsWith(r,'RS')
                r='RSC';
            end
            if startsWith(r,'VISpm')
                r='VISm';
            end
            try
                
                
                dat = nanmean(tmp_factors(:,iidx))';
                dat(~enough_data_idx(iidx))=nan;
                dat2 = nanmean(tmp_factors2(:,iidx))';
                dat2(~enough_data_idx2(iidx))=nan;
                stab = nanmean(a.all_stability(:,iidx))';
                
                tmp = [stab,dat,ones(size(dat))*iF,ones(size(dat))*aid,dat2]; %stability, factors,recording_id, firing rate
                tmpgc = [sum(~isnan(tmp_factors(:,iidx)),1)', ones(nnz(iidx),1)*size(a.all_stability,1)];
                if ismember(r,fieldnames(STABILITY))
                    STABILITY.(r) = cat(1,STABILITY.(r),tmp);
                    GoodCells.(r) = cat(1,GoodCells.(r),tmpgc);
                else
                    STABILITY.(r)=tmp;
                    GoodCells.(r)=tmpgc;
                end
            catch ME
                disp(ME.message)
            end
        end
    catch ME
        disp(ME.message)
    end
    
end