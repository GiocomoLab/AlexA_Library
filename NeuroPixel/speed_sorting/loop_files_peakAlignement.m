
ops.BinWidth = 1;

ops.edges = 0:ops.BinWidth:400;
ops.search_range=22/ops.BinWidth:1:378/ops.BinWidth; % in bins

ops.midpoints = ops.edges(1:end-1)+(ops.edges(1)+ops.edges(2))/2;
ops.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = [];
ops.n_preceeding = 10;
ops.TimeBin = 0.02;
ops.speedWindow = [-10 -1]; % in cm

fi = gausswin(11);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = false;
ops.idx = ops.search_range;% in bins, for calculating corr

OAK='/oak/stanford/groups/giocomo/';

ops_shifts.factors = -1:0.01:1;
ops_shifts.edges = 0:1:400;
ops_shifts.nBins = numel(ops_shifts.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops_shifts.trials = 3:20;
ops_shifts.TimeBin = 0.02;
ops_shifts.idx = [10:1:390]/1;% in bins
fi = gausswin(5);
fi=fi'/sum(fi);
ops_shifts.filter = fi;
ops_shifts.plotfig = false;
ops.ops_shifts = ops_shifts;


%%
% data_dir=fullfile('F:','NP_DATA');
% %session_name = {'AA5_190809_gain_1'};
% filenames = {};
% sn = dir(fullfile(data_dir,'*.mat'));
% for iS = 1:numel(sn)
%     if ~(contains(sn(iS).name,'mismatch') || contains(sn(iS).name,'playback') || contains(sn(iS).name,'dark'))
%         filenames{end+1}=sn(iS).name(1:end-4);
%     end
% end
gain = 0.8;
contrast = 100;
region = 'MEC';
[filenames,triggers] = getFilesCriteria(region,contrast,gain,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
savepath='/oak/stanford/groups/giocomo/attialex/images_peakLong_gausswin_5';
if ~isfolder(savepath)
    mkdir(savepath)
end
%%
if isfile(fullfile(savepath,[region '.mat']))
    resp = input([region 'output file already exists,continue [1/0]']);
    if resp==0
        disp('abort, load ouput file')
        load(fullfile(savepath,[region '.mat']))
        return
    end
    
end
save(fullfile(savepath,'parameters.mat'),'ops')
%%

p=gcp('nocreate');
if isempty(p)
    parpool(12);
end
%%
x_vec =-20/ops.BinWidth:1:20/ops.BinWidth;

%%
plotFigures=false;
output = cell(numel(filenames),1);
parfor iF=1:numel(filenames)
    try
        %data = load(fullfile(data_dir,filenames{iF}));
        data = load(filenames{iF});
        ops_temp = ops;
        data_out=cell(2,1);
        for iRep = 1:min(2,numel(triggers{iF}))
            preceeding_trials = triggers{iF}(iRep)+[-ops_temp.n_preceeding:-1];
            if ~all(data.trial_gain(preceeding_trials)==1 & data.trial_contrast(preceeding_trials)==contrast)
                continue
            end
            ops_temp.trials=triggers{iF}(iRep)+[-ops_temp.n_preceeding:3];
            data_out{iRep} = findPeakAlignement(data,ops_temp);
            data_out{iRep}.trials = ops_temp.trials;
            test_trials = triggers{iF}(iRep)+[4:13];
            ops_temp.ops_shifts.trials=test_trials;
            [correlation_shifted,correlation_noshift]=test_alignement(data,ops_temp.ops_shifts,data_out{iRep}.factors,data_out{iRep}.CID);
            data_out{iRep}.correlation_shifted = correlation_shifted;
            data_out{iRep}.correlation_noshift = correlation_noshift;
            if plotFigures
                fig = figure('Position',[680   758   560   220],'visible','off','Renderer','Painters');
                
                for iC=1:numel(data_out{iRep}.region)
                    tmp_d = data_out{iRep}.correlation_shifted(iC)-data_out{iRep}.correlation_noshift(iC);

                    if isnan(data_out{iRep}.max_ind(iC)) || ~(data_out{iRep}.stability(iC)>.2) || ~startsWith(data_out{iRep}.region{iC},region)
                        continue
                    end
                    trial_speed = getSpeedAroundPoint(data_out{iRep}.speed(:,1),data.posx,data.trial,ops_temp,data_out{iRep}.max_ind(iC),ops_temp.speedWindow);
                    trial_speed_gain = trial_speed(end-3:end);
                    
                    trial_speed = trial_speed(1:end-4);
                    [~,sidx]=sort(trial_speed,'descend');
                    [~,sidx_gain]=sort(trial_speed_gain,'descend');
                    tmp_rank = 1:ops_temp.n_preceeding;
                    tmp_rank(sidx)=tmp_rank;
                    tmp_rank_gain = ops_temp.n_preceeding+[1:4];
                    tmp_rank_gain(sidx_gain)=tmp_rank_gain;
                    tmp_rank = [tmp_rank, tmp_rank_gain];
                    tmp = data_out{iRep}.allSpikes{iC}(:,2);
                    for iT=1:numel(tmp)
                        tmp(iT)=tmp_rank(tmp(iT));
                    end
                    
                    
                    subplot(2,3,1)
                    scatter(data_out{iRep}.allSpikes{iC}(:,1),tmp,15,get_color(1,contrast),'.')
                    gain_idx = tmp>ops.n_preceeding;
                    hold on
                    scatter(data_out{iRep}.allSpikes{iC}(gain_idx,1),tmp(gain_idx),15,get_color(gain,contrast),'.')
                    xlim(data_out{iRep}.max_ind(iC)+[-60 40] )
                    ylim([0 23])
                    xline(data_out{iRep}.max_ind(iC));
                    box off
                    title(data_out{iRep}.region{iC})
                    
                    subplot(2,3,2)
                    scatter(trial_speed(sidx),[1:ops_temp.n_preceeding],15,get_color(1,contrast),'.')
                    
                    hold on
                    scatter(trial_speed_gain(sidx_gain),ops_temp.n_preceeding+[1:4],15,get_color(gain,contrast),'.')
                    box off
                    xlabel('cm/s')
                    %plot the corrected spikes, they are in column 3
                    subplot(2,3,4)
                    scatter(data_out{iRep}.allSpikes{iC}(:,3),tmp,15,get_color(1,contrast),'.')
                    hold on
                    scatter(data_out{iRep}.allSpikes{iC}(gain_idx,3),tmp(gain_idx),15,get_color(gain,contrast),'.')
                    xlabel(sprintf('%.2f, diff: %.2f',data_out{iRep}.factors(iC),tmp_d))
                    xlim(data_out{iRep}.max_ind(iC)+[-60 40] )
                    ylim([0 23])
                    xline(data_out{iRep}.max_ind(iC));
                    if ~isempty(data_out{iRep}.subregion)
                        title(data_out{iRep}.subregion{iC})
                    end
                    
                    box off
                    subplot(2,3,[3 6])
                    plot(x_vec,data_out{iRep}.all_fast(iC,:))
                    hold on
                    plot(x_vec,data_out{iRep}.all_slow(iC,:))
                    plot(x_vec,data_out{iRep}.all_gain(iC,:));
                    legend({'fast','slow','gain'});
                    grid on
                    box off
                    [~,sn,~]=fileparts(filenames{iF});
                    im_save_path = fullfile(savepath,sprintf('%s_%.1f',sn,gain));
                    if ~isfolder(im_save_path)
                        mkdir(im_save_path)
                    end
                    impath = fullfile(im_save_path,sprintf('%d_%d.png',data_out{iRep}.CID(iC),iRep));
                    saveas(gcf,impath);
                    impath = fullfile(im_save_path,sprintf('%d_%d.pdf',data_out{iRep}.CID(iC),iRep));
                    saveas(gcf,impath);
                    %pause
                    clf
                end
                close(fig)
            end
        end
        output{iF}=data_out;
    catch ME
        disp(ME.message)
        
        disp('whoopsie')
    end
end
%%
save(fullfile(savepath,region),'output')

%%
