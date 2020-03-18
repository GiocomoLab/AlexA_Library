%% struct_names
dataset = 'filtered_gaincontrast_space11_speed5';
path = ['Z:\giocomo\attialex\speed_' dataset];
filenames = dir(fullfile(path,'*.mat'));

%% 

FACTOR_HIGH = struct();
FACTOR_LOW = struct();
STABILITY = struct();
DEPTH = struct();
SITE = struct();
for iF=1:numel(filenames)
    try
    filename=filenames(iF).name;
    a=load(fullfile(path,filename));
    [unique_regions,~,ridx]=unique(a.region);
    for iR=1:numel(unique_regions)
        
        r=unique_regions{iR};
        iidx = ridx==iR;
        if startsWith(r,'RS')
            r='RSC';
        end
        if startsWith(r,'VISpm')
            r='VISm';
        end
        try
            tmp_high = nanmean(a.all_factors(a.chunk_contrast==100,:));
            tmp_stab = nanmean(a.all_stability);
            tmp_low = nanmean(a.all_factors(a.chunk_contrast==10,:));
        if ismember(r,fieldnames(FACTOR_HIGH))
        FACTOR_HIGH.(r) = cat(2,FACTOR_HIGH.(r),[tmp_high(iidx)]);
        DEPTH.(r) = cat(2,DEPTH.(r),a.depth(iidx));
        FACTOR_LOW.(r) = cat(2,FACTOR_LOW.(r),[tmp_low(iidx)]);
        STABILITY.(r) = cat(2,STABILITY.(r),tmp_stab(iidx));
        SITE.(r) = cat(1,SITE.(r),iF*ones(nnz(iidx),1));
        else
            FACTOR_HIGH.(r)=[tmp_high(iidx)];
            FACTOR_LOW.(r) = tmp_low(iidx);
            STABILITY.(r) = tmp_stab(iidx);
            DEPTH.(r) = a.depth(iidx);
            SITE.(r) = iF*ones(nnz(iidx),1);
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
ops = load(fullfile(path,'parameters'));
ops = ops.ops;
fn = fieldnames(FACTOR_HIGH);
fn = {'VISp','MEC'};
histfig=figure;
for iF=1:numel(fn)
    dat = FACTOR_HIGH.(fn{iF});
    %figure
    if ismember(fn{iF},{'MEC','ECT'})
        mult = -1;
    else
        mult = 1;
    end
    %idx = ma>0.1; %correlaton at 0 lag
    
    tmp_fact = dat;
    tmp_stab = STABILITY.(fn{iF});
    idx = tmp_stab>0.2;
    subplot(5,3,[1 4 7 10]+iF-1)
    plot(tmp_fact(idx),DEPTH.(fn{iF})(idx)*mult,'.')
    set(gca,'YDir','reverse')
    title(fn{iF})
    xlim([-.3 .3])
    xlabel('speed correction lag [s]')
    grid on
    subplot(5,3,12+iF)
    histogram(tmp_fact(idx),'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5],'EdgeColor','none')
    xlim([-.3 .3])
    xlabel('speed correction lag [s]')
end
    
   
%%
ops = load(fullfile(path,'parameters'));
ops = ops.ops;
fn = fieldnames(FACTOR_HIGH);
fn = {'VISp','MEC'};
histfig=figure;
for iF=1:numel(fn)
    dat = FACTOR_HIGH.(fn{iF});
    %figure
    
    %idx = ma>0.1; %correlaton at 0 lag
    
    tmp_fact = dat;
    tmp_stab = STABILITY.(fn{iF});
    idx = tmp_stab>.2;
    figure(histfig)
    subplot(numel(fn),2,iF*2-1)
    hold on
    histogram(tmp_fact(idx),'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])
    histogram(FACTOR_LOW.(fn{iF})(idx),'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])
    subplot(numel(fn),2,iF*2)
    hold on
    [f,xi]=ksdensity(tmp_fact(idx));
    [f_low,xi_low]=ksdensity(FACTOR_LOW.(fn{iF})(idx));
    plot(xi,f)
    hold on
    plot(xi_low,f_low)
    legend({'100','10'})
end


    xlabel('speed correction lag [s]')
    
%%
ops = load(fullfile(path,'parameters'));
ops = ops.ops;
fn = fieldnames(FACTOR_HIGH);
fn = {'VISp','MEC'};
histfig=figure;
lims=[-.4 .6];
for iF=1:numel(fn)
    subplot(1,numel(fn),iF)
    idx = STABILITY.(fn{iF})>.2;
    scatter(FACTOR_HIGH.(fn{iF})(idx),FACTOR_LOW.(fn{iF})(idx))
    axis image
    xlim(lims)
    ylim(lims)
    grid on
    xlabel('shift high contrast')
    ylabel('shift low contrast')
end
%%
ops = load(fullfile(path,'parameters'));
ops = ops.ops;
fn = fieldnames(FACTOR_HIGH);
fn = {'VISp','MEC'};
histfig=figure;
lims=[-.4 .6];
for iF=1:numel(fn)
    subplot(1,numel(fn),iF)
    idx = STABILITY.(fn{iF})>.2;
    [a,b,c]=unique(SITE.(fn{iF}));
    tmp_high=nan(numel(a),1);
    tmp_low = tmp_high;
    for iS=1:numel(a)
        idx_s = c==iS & idx';
        tmp_high(iS)=nanmean(FACTOR_HIGH.(fn{iF})(idx_s));
        tmp_low(iS) = nanmean(FACTOR_LOW.(fn{iF})(idx_s));
    end
    scatter(tmp_high,tmp_low)
    axis image
    xlim(lims)
    ylim(lims)
    grid on
    xlabel('shift high contrast')
    ylabel('shift low contrast')
end
