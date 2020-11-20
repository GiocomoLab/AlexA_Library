%% struct_names
%first run loop_findshifts to generate the data
path = '/Volumes/Samsung_T5/attialex/speed_filtered_correctedData';
path = '/Volumes/Crucial X8/attialex/speed_filtered_correctedData_shortidx2';

fig_path='/Volumes/Samsung_T5';
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
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile('/Volumes/Crucial X8/attialex/NP_DATA_corrected'));
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
                
                enough_data_idx = sum(~isnan(tmp_factors),1)>2; %only use data where there are more than 2 stable blocks
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
                
                stab = nanmean(a.all_stability(:,iidx))';
                
                tmp = [stab,dat,ones(size(dat))*iF,ones(size(dat))*aid]; %stability, factors,recording_id, firing rate
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

%%
reg = {'MEC','VISp','RSC'};
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
AV={};
AV_ALL={};
X=[];
G=[];
AID=[];
figure

for iR = 1:numel(reg)
    Xt=[];
    Gt=[];
    AID=[];
    counter=0;
    subplot(2,3,iR)
    hold on
    data_this = STABILITY.(reg{iR});
    %data_this(data_this<=-.15)=nan;
    [a,b]=unique(data_this(:,3));
    averages  = zeros(size(a));
    for iP = 1:numel(a)
        idx = data_this(:,3)==a(iP) & ~isnan(data_this(:,2));
        averages(iP)=nanmedian((data_this(idx,2)))*-1;
        X(end+1)=averages(iP);
        G(end+1)=iR;
        if nnz(idx)<5;%nnz(idx)<nnz(data_this(:,3)==a(iP))*.1
            averages(iP)=nan;
            
        else
            counter=counter+1;
            p=signrank(data_this(idx,2)*-1);
            if p<0.05
                filled = 'filled'
            else
                filled=[];
            end
            scatter(data_this(idx,2)*-1,ones(1,nnz(idx))*averages(iP),filled)
            Xt=cat(1,Xt,data_this(idx,2)*-1);
            Gt=cat(1,Gt,ones(nnz(idx),1)*counter);
            AID = cat(1,AID,data_this(idx,4));
        end
    end
    [a,b]=sort(averages(~isnan(averages)));
    ranking = 1:numel(b);
    ranking(b)=ranking;
    ylim([-.1 0.3])
    subplot(2,3,iR+3)
    Gt_sorted = ranking(Gt);
    AID_sorted = AID;
    for ii=1:numel(b)
        take_idx = Gt==ii;
        put_idx  = Gt_sorted==ranking(ii);
        AID_sorted(put_idx)=AID(take_idx);
    end
        boxplot(Xt,Gt_sorted,'PlotStyle','compact','ColorGroup',AID_sorted)
    ylim([-.3 .3])
    title(regions{iR})
    ylabel('delay')
end

