fi = dir('/Volumes/T7/attialex/speed_filtered_Playback/AA*back*Towers*.mat');


%%
gain_names={};
pb_names = {};
for iF = 1:numel(fi)
    
    data=load(fullfile(fi(iF).folder,fi(iF).name),'anatomy');
    if numel(fields(data)) ==0
        continue
    end
    idx = strfind(fi(iF).name,'_');
    newStr = [fi(iF).name(1:idx(2)), '*.mat'];
    fa = dir(fullfile(fi(iF).folder,newStr))
    pp = {fa(:).name}
    pot = ~contains(pp,{'mismatch','play','dark'})
    nPot = nnz(pot)
    if nPot>0
        idx = find(nPot,1);
        gain_names{end+1}=pp{idx};
        pb_names{end+1}=fi(iF).name;
    end
    
    
end

%%
opt = load_default_opt;
NBlocks=[];
Frac_Good = [];
STABILITY = struct();
ANIMALS = struct();
animal_counter = 0;
for iF=1:numel(gain_names)
    data = load(fullfile('/Volumes/T7/attialex/speed_filtered_correctedData_shortidx2',gain_names{iF}));
    data_pb = load(fullfile('/Volumes/T7/attialex/speed_filtered_Playback',pb_names{iF}));
    [unique_regions,~,ridx]=unique(data.region);

        
        parts=strsplit(pb_names{iF},'_');
        animal = parts{1};
        if ~isfield(ANIMALS,animal)
            animal_counter = animal_counter +1;
            ANIMALS.(animal) = animal_counter;
        end
        aid = ANIMALS.(animal);
    
    
    NBlocks = cat(1,NBlocks,size(data.all_stability,1));
    stab = data.all_stability;
    stable_blocks = stab>=opt.stab_thresh;
    tmp_factors = data.all_factors;
    tmp_factors(~stable_blocks)=nan; %set shifts to nan where stability <.5
    stab_pb = data_pb.all_stability;
    stable_blocks_pb = stab_pb>-opt.stab_thresh;
    tmp_factors_pb = data_pb.all_factors;
    tmp_factors_pb(~stable_blocks_pb)=nan;

    enough_data_idx = sum(~isnan(tmp_factors),1)>2; %only use data where there are more than 2 stable blocks
    tmp_idx = startsWith(data.region,regions);
    Frac_Good = cat(1,Frac_Good,[sum(~isnan(tmp_factors(:,tmp_idx)),1)', ones(nnz(tmp_idx),1)*size(data.all_stability,1)]);
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
            
            
            dat = nanmean(tmp_factors(:,iidx),1)';
            dat(~enough_data_idx(iidx))=nan;
            
            stab = nanmean(data.all_stability(:,iidx),1)';
            
            tmp = [stab,dat,ones(size(dat))*iF,ones(size(dat))*aid,nanmean(tmp_factors_pb(:,iidx),1)']; %stability, factors,recording_id, firing rate
            tmpgc = [sum(~isnan(tmp_factors(:,iidx)),1)', ones(nnz(iidx),1)*size(data.all_stability,1)];
            if ismember(r,fieldnames(STABILITY))
                STABILITY.(r) = cat(1,STABILITY.(r),tmp);
                GoodCells.(r) = cat(1,GoodCells.(r),tmpgc);
            else
                STABILITY.(r)=tmp;
                GoodCells.(r)=tmpgc;
            end
        catch ME
            disp(ME.message)
            if startsWith(ME.message,'Dim')
                keyboard
            end
        end
    end
    
end
%%
figure
subplot(1,2,1)
scatter(STABILITY.VISp(:,2),STABILITY.VISp(:,end))
refline(1,0)
axis image
xlabel('Delay Gain Session')
ylabel('Delay PB Session')
title('VISp')
subplot(1,2,2)
scatter(STABILITY.RSC(:,2),STABILITY.RSC(:,end))
refline(1,0)
axis image
xlabel('Delay Gain Session')
ylabel('Delay PB Session')
title('RSC')

