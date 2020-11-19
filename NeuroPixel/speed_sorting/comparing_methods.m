%% struct_names
path = '/Volumes/Crucial X8/attialex/speed_filtered_correctedData_shortidx2';
path1 = '/Users/attialex/speed_filtered_correctedData';
path1='/Users/attialex/speed_filtered_correctedData_oldOcc';
%path =  '/Volumes/Samsung_T5/tbtxcorr_08';
%path = '/Users/attialex/speed_filtered_correctedData';
fig_path='/Volumes/Samsung_T5';
%fig_path ='C:\figs\campbell_attinger\fig4_sensory_delay';
filenames = dir(fullfile(path,'*.mat'));
ops = load(fullfile(path,'parameters.mat'));
ops = ops.ops;
opt = load_default_opt;
regions = {'MEC','VISp','RSC'};

%%

STABILITY = struct();
STAB_FL = struct();
GoodCells = struct();
NBlocks=[];
Frac_Good=[];
for iF=1:numel(filenames)
    try
        filename=filenames(iF).name;
        
        a=load(fullfile(path,filename));
        b=load(fullfile(path1,filename));
        [unique_regions,~,ridx]=unique(a.region);
        NBlocks = cat(1,NBlocks,size(a.all_stability,1));
        stab = a.all_stability;
        stable_blocks = stab>=opt.stab_thresh;
        tmp_factors = a.all_factors;
        tmp_factors(~stable_blocks)=nan; %set shifts to nan where stability <.5
        %tmp_factors_mgc = b.all_factors_mgc;
        tmp_factors_mgc = b.all_factors;
        tmp_factors_mgc(~stable_blocks)=nan;
        enough_data_idx = sum(~isnan(tmp_factors),1)>=2; %only use data where there are more than 2 stable blocks
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
                
                dat_fl = [a.all_factors(1,iidx)',a.all_factors(end,iidx)'];
                stab_fl = [a.all_stability(1,iidx)',a.all_stability(end,iidx)'];
                tmp_fl = cat(2,dat_fl,stab_fl);
                dat = nanmean(tmp_factors(:,iidx))';
                dat(~enough_data_idx(iidx))=nan;
                dat_mgc = nanmean(tmp_factors_mgc(:,iidx))';
                dat_mgc(~enough_data_idx(iidx)) = nan;
                
                stab = nanmean(a.all_stability(:,iidx))';
                tmp = [stab,dat,ones(size(dat))*iF,dat_mgc]; %stability, factors,recording_id, firing rate
                tmpgc = [sum(~isnan(tmp_factors(:,iidx)),1)', ones(nnz(iidx),1)*size(a.all_stability,1)];
                if ismember(r,fieldnames(STABILITY))
                    STABILITY.(r) = cat(1,STABILITY.(r),tmp);
                    GoodCells.(r) = cat(1,GoodCells.(r),tmpgc);
                    STAB_FL.(r)=cat(1,STAB_FL.(r),tmp_fl);
                else
                    STABILITY.(r)=tmp;
                    GoodCells.(r)=tmpgc;
                    STAB_FL.(r)=tmp_fl;
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
figure('Position',[440   464   625   334])
regions = {'MEC','VISp','RSC'};
for ii=1:3
    subplot(2,3,ii)
    scatter(STABILITY.(regions{ii})(:,2)*-1,STABILITY.(regions{ii})(:,4)*-1)
    xlabel('custom')
    ylabel('corrMat')
    axis image
    xlim([-.3 .3])
    ylim([-.3 .3])
    refline(1,0)
    title(regions{ii})
    subplot(2,3,ii+3)
    histogram(STABILITY.(regions{ii})(:,2)*-1)
    hold on
    histogram(STABILITY.(regions{ii})(:,4)*-1)
end
legend('Custom','CorrMat')

%%
reg = {'VISp','RSC','MEC','RHP','ECT'};
reg = {'MEC','VISp','RSC'}
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
AV={};
AV_ALL={};
X=[];
G=[];
figure
for iR = 1:numel(reg)
    data_this = STABILITY.(reg{iR});
    %data_this(data_this<=-.15)=nan;
    [a,b]=unique(data_this(:,3));
    averages  = zeros(numel(a),2);
    for iP = 1:numel(a)
        idx = data_this(:,3)==a(iP);
        averages(iP,:)=nanmean((data_this(idx,[2,4]))*-1);
        
        X(end+1)=averages(iP);
        G(end+1)=iR;
        if nnz(idx)<nnz(data_this(:,3)==a(iP))*.1
            averages(iP)=nan;
        end
    end
    AV_ALL.(reg{iR}) = averages;
    AV{iR}=averages;
    
end
%%
figure('Position',[440   218   628   189]);

for ii=1:numel(regions)
    subplot(1,3,ii)
    scatter(AV_ALL.(regions{ii})(:,1),AV_ALL.(regions{ii})(:,2))
    axis image
    xlim([-.1 .35])
    ylim([-.1 .35])
    refline(1,0)
    xlabel('custom')
    ylabel('corrMat')
    title(regions{ii})
    
end
%%
figure
for ii=1:3
    subplot(1,3,ii)
    val_idx = all(STAB_FL.(regions{ii})(:,3:4)>.5,2);
    scatter(STAB_FL.(regions{ii})(val_idx,1)*-1,STAB_FL.(regions{ii})(val_idx,2)*-1)
    axis image
    xlim([-.3 .3])
    ylim([-.3 .3])
    refline(1,0)
    title(regions{ii})
    disp(signrank(STAB_FL.(regions{ii})(val_idx,1)*-1,STAB_FL.(regions{ii})(val_idx,2)*-1))
end
