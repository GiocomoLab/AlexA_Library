ops = load_default_opt;
ops.trial_range = [-6:9];
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.max_lag = 30;
ops.maxLag = ops.max_lag;
OAK='/oak/stanford/groups/giocomo/';
%OAK = '/Volumes/Samsung_T5';
gains = [0.5, 0.6, 0.7, 0.8];
%%
gain = 0.6;
contrast = 100;
regions = {'MEC'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(OAK,'attialex','NP_DATA_corrected'));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end
savepath = fullfile(OAK,'attialex','tbtxcorr_range_2');
fig_save_path = fullfile(savepath,'images');

shiftDir = fullfile(OAK,'attialex','speed_filtered_correctedData');
if ~isfolder(savepath)
    mkdir(savepath)
end
if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end
shift_ops = load(fullfile(shiftDir,'parameters.mat'));
shift_ops = shift_ops.ops;
save(fullfile(OAK,'attialex','parameters.mat'),'ops');
%%
p = gcp('nocreate');
if isempty(p)
    p = parpool(12);
end
%%

%%
parfor iF=1:numel(filenames)
   fig=figure('Renderer','Painters','Position',[440   611   877   187],'Visible','off');
cmap = cbrewer('seq','BuPu',20);

colormap(cmap) 
    try
        [~,sn]=fileparts(filenames{iF});
        mkdir(fullfile(fig_save_path,sn))
        
        
        data = load(filenames{iF});
        
        if isfield(data.anatomy,'parent_shifted')
            reg = data.anatomy.parent_shifted;
        else
            reg = data.anatomy.cluster_parent;
        end
        if iscolumn(reg)
            reg = reg';
        end
        
        reg=reg(data.sp.cgs==2);
        reg_orig = data.anatomy.cluster_parent((data.sp.cgs==2));
        if iscolumn(reg_orig)
            reg_orig = reg_orig';
        end
        
        
        
        %calculate trial by trial correlation across all bins
        trials=1:max(data.trial);
        ops_here = ops;
        ops_here.trials = trials;
        cellID = data.sp.cids(data.sp.cgs==2);
        
        [corrMat,frMat,shiftMat]=trialCorrMat(cellID,trials,data,ops);
        ons05=strfind(data.trial_gain'==0.5,[0 1])
        if numel(ons05)<2
            continue
        end
        for iC=1:size(frMat,1)
            good=0;
            for iG=1%1:numel(gains)
                onsets = strfind(data.trial_gain' == gains(iG),[0 1])+1;
                pre_map = zeros(2,200);
                pre_map_all = zeros(2,10,200);
                
                gain_map = pre_map;
                pre_stab = zeros(1,2);
                gain_stab = zeros(1,2);
                
                for iO=1:2
                    onset = onsets(iO);
                    pre_idx = onset-6:onset-1;
                    all_idx = onset-6:onset+3;
                    gain_idx = onset+(0:3);
                    pre_map_all(iO,:,:)=squeeze(frMat(iC,all_idx,:));
                    pre_map(iO,:) = squeeze(mean(frMat(iC,pre_idx,:),2));
                    gain_map(iO,:)=squeeze(mean(frMat(iC,gain_idx,:),2));
                    pre_stab(iO) = nanmean(nanmean(corrMat(iC,pre_idx,pre_idx),2),3);
                    gain_stab(iO) = nanmean(nanmean(corrMat(iC,gain_idx,gain_idx),2),3);
                end
                
                if all(pre_stab>.5) && all(gain_stab>.5)
                    good=good+1;
                    subplot(3,4,iG)
                    imagesc([pre_map(1,:);gain_map(1,:);pre_map(2,:);gain_map(2,:)])
                    
                    subplot(3,4,iG+4)
                       trials = onsets(1)+(-6:3);

                    imagesc(ops.xbincent,trials,squeeze(pre_map_all(1,:,:)))
                    hold on
                    for tr = trials
                        patch([-15 0 0 -15],[tr tr tr+1 tr+1]-0.5,get_color(data.trial_gain(tr),100),...
                            'EdgeColor',get_color(data.trial_gain(tr),100));
                    end
                    xlim([-15 400]);
                    
                    subplot(3,4,iG+8)
                                        trials = onsets(2)+(-6:3);

                    imagesc(ops.xbincent,trials,squeeze(pre_map_all(2,:,:)))
                    
                    for tr = trials
                        patch([-15 0 0 -15],[tr tr tr+1 tr+1]-0.5,get_color(data.trial_gain(tr),100),...
                            'EdgeColor',get_color(data.trial_gain(tr),100));
                    end
                    xlim([-15 400]);
                    
                end
                
            end
            %colormap(cmap)
            
            %pause
            if good==1
                cluID = cellID(iC);
                saveas(fig,fullfile(fig_save_path,sn,sprintf('%d.png',cluID)))
                saveas(fig,fullfile(fig_save_path,sn,sprintf('%d.pdf',cluID)))
            end
            clf
            
        end
        
        close(fig)
        
        
        
        data_out = matfile(fullfile(savepath,sprintf('%s',sn)),'Writable',true);
        data_out.corrMat = corrMat;
        data_out.shiftMat = shiftMat;
        
        
        
        data_out.region = reg;
        
        
        
        data_out.trials = trials;
        data_out.gain = data.trial_gain;
        
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end

