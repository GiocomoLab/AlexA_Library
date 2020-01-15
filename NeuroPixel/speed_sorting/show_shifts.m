%% struct_names
dataset = 'shiftNoFilter';
path = ['Z:\giocomo\attialex\speed_' dataset];
filenames = dir(fullfile(path,'*.mat'));

%% 

STABILITY = struct();

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
        if startsWith(r,'VIS')
            r='VISp';
        end
        try
        if ismember(r,fieldnames(STABILITY))
        STABILITY.(r) = cat(1,STABILITY.(r),a.stability(iidx,:));
        else
            STABILITY.(r)=a.stability(ridx,:);
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
fn = fieldnames(STABILITY);
fn = {'MEC','VISp'};
histfig=figure;
for iF=1:numel(fn)
    dat = STABILITY.(fn{iF});
    figure
    [ma,mi]=max(dat,[],2);
    dif = ma-dat(:,26);
    idx = ma>0.1;
    plot(ops.factors(mi),ma,'.')
    figure(histfig)
    hold on
    histogram(ops.factors(mi(idx)),'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])

end
    legend(fn)
    xlabel('speed correction lag [s]')
    
%%
fn = fieldnames(STABILITY);

N = zeros(numel(fn));
for iF=1:numel(fn)
    N(iF)=size(STABILITY.(fn{iF}),1);
end

[a,sid]=sort(N,'descend');
figure
for iF=1:5
    
    dat = STABILITY.(fn{sid(iF)});
    subplot(5,2,iF*2)
    
    [ma,mi]=max(dat,[],2);
    idx = ma>0.01;
    plot(ops.factors(mi),ma,'.')
    subplot(5,2,iF*2-1)
    hold on
    histogram(ops.factors(mi(idx)),'Normalization','probability','BinEdges',[ops.factors-mean(diff(ops.factors))*.5])
        legend(fn{sid(iF)})
    xlabel('speed correction lag [s]')
end

    
    