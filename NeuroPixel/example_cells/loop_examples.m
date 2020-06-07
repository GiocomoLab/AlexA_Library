% firing rate maps for example cells in small gain manip

paths = struct;
paths.data = '/Volumes/Samsung_T5/attialex/NP_Data_corrected';
paths.figs = '/Volumes/Samsung_T5/attialex/images/example_cells4'; % where to save figs
if ~isfolder(paths.figs)
    mkdir(paths.figs)
end
opt = load_default_opt;
opt.num_tr_bl = 6;

%%
gain = 0.5;
contrast = 100;
regions = {'MEC'};
filenames = {};
triggers = {};
OAK = paths.data;
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
paths.figs = fullfile(paths.figs,regions{1});
if ~isfolder(paths.figs)
    mkdir(paths.figs)
end

%%


for iF = 1:numel(filenames)
    dat = load(filenames{iF});

    for iT=1:numel(triggers{iF})
        trials = triggers{iF}(iT)-6+(0:15);

    
    if isfield(dat,'anatomy')
        if isfield(dat.anatomy,'cluster_parent')
            region = dat.anatomy.cluster_parent;
        else
            region = dat.anatomy.parent_shifted;
        end
    else
        continue
    end
    if iscolumn(region)
        region = region';
    end

    good_idx = dat.sp.cgs==2 & startsWith(region,regions{1});
    good_cells = dat.sp.cids(good_idx)';
    STAB = zeros(size(good_cells));
    %FRMAT = zeros(16,200,numel(good_cells));
    region = region(good_idx);
    if isfield(dat.anatomy,'depth')
    depth  = dat.anatomy.depth(good_idx);
    else
        depth=zeros(size(STAB));
    end
    [corrMat,frMatAll]=trialCorrMat(good_cells,trials,dat,opt);
    STAB=nanmean(nanmean(corrMat(:,1:6,1:6),2),3);

    % compute firing rate and firing rate correlation matrices
    [~,sid]=sort(STAB,'descend','MissingPlacement','Last');
    [~,sn] = fileparts(filenames{iF});
    % make fig
    n2plot = min(10,size(frMatAll,3));
   hfig=figure('Renderer','Painters','Position',[440   609   452   189],'Visible','off');
cmap = cbrewer('seq','BuPu',20);

colormap(cmap)
    for iC=1:n2plot
        frMat = squeeze(frMatAll(sid(iC),:,:));
        hold on
        cluID = good_cells(sid(iC));
        rr=region{sid(iC)};
        dd = depth(sid(iC));
        celltitle_save = sprintf('%s_c%d_%s_d=%d_tr%d-%d_stab=%0.2f',sn,cluID,...
            rr,dd,...
            trials(1),trials(end),STAB(sid(iC)));
        celltitle_fig = sprintf('%s c%d\n%s d=%d stab_bl=%0.2f',sn,cluID,...
            rr,dd,STAB(sid(iC)));
        
%         hfig(iC) = figure('Position',[200 200 400 400],'Visible','off'); hold on;
%         hfig(iC).Name = celltitle_save;

            imagesc(opt.xbincent,trials,squeeze(frMat));
    text(200,max(trials)+1.5,sprintf('Stab: %.2f, Peak FR; %d Hz',STAB(sid(iC)),round(max(frMat,[],'all'))))
    %title(celltitle_fig,'Interpreter','none');
    ylim([min(trials)-.5 max(trials)+.5]);
    yticks([min(trials) max(trials)]);
    ylabel('Trial');
    xticks([1 200 400]);
    xticklabels({'0', 'cm','400'});
%     %xlabel('cm');
%         title(celltitle_fig,'Interpreter','none');
%         ylim([min(trials) max(trials)]);
%         yticks([min(trials) max(trials)]);
%         ylabel('trial');
%         xlim([0 400]);
%         xticks([0 200 400]);
%         xticklabels([0 200 400]);
%         xlabel('cm');
%         
        %patches indicating gain value
        for tr = trials
            patch([-15 0 0 -15],[tr tr tr+1 tr+1]-0.5,get_color(dat.trial_gain(tr),100),...
                'EdgeColor',get_color(dat.trial_gain(tr),100));
        end
        for tr = 0.5+[trials(opt.num_tr_bl) trials(opt.num_tr_bl+opt.num_tr_gc)]
            plot([0 400],[tr tr],'w-');
        end
        xlim([-15 400]);
        
        
        % save figs
        
        
    
    [~,sn]=fileparts(filenames{iF});
    saveas(gcf,fullfile(paths.figs,sprintf('%s_%d_%d.pdf',sn,cluID,iT)));
    clf
    end
        close(hfig)
    end

end