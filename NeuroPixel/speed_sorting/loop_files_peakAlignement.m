
ops.factors = -.25:0.01:.25;
ops.BinWidth = 1;

ops.edges = 0:ops.BinWidth:400;
ops.search_range=22/ops.BinWidth:1:378/ops.BinWidth; % in bins

ops.midpoints = ops.edges(1:end-1)+(ops.edges(1)+ops.edges(2))/2;
ops.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = 3:24;
ops.TimeBin = 0.02;
ops.speedWindow = [-10 -1]; % in cm

fi = gausswin(11);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = false;
ops.idx = ops.search_range;% in bins, for calculating corr

OAK='/oak/stanford/groups/giocomo/';

ops_shifts.factors = -.25:0.01:.25;
ops_shifts.edges = 0:2:400;
ops_shifts.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops_shifts.trials = 3:20;
ops_shifts.TimeBin = 0.02;
ops_shifts.idx = [10:2:390]/2;% in bins
fi = gausswin(5);
fi=fi'/sum(fi);
ops_shifts.filter = fi;
ops_shifts.plotfig = true;
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
region = 'ECT';
[filenames,triggers] = getFilesCriteria(region,contrast,gain,'/oak/stanford/groups/giocomo/attialex/NP_DATA');
savepath='/oak/stanford/groups/giocomo/attialex/images_peak';
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
            
            ops_temp.trials=triggers{iF}(iRep)+[-18:3];
            data_out{iRep} = findPeakAlignement(data,ops_temp);
            
            if plotFigures
                fig = figure('Position',[680   758   560   220],'visible','off','Renderer','Painters');
                
                for iC=1:numel(data_out{iRep}.region)
                    if isnan(data_out{iRep}.max_ind(iC)) || ~(data_out{iRep}.stability(iC)>.2) || ~startsWith(data_out{iRep}.region{iC},region)
                        continue
                    end
                    trial_speed = getSpeedAroundPoint(data_out{iRep}.speed(:,1),data.posx,data.trial,ops,data_out{iRep}.max_ind(iC),ops.speedWindow);
                    trial_speed_gain = trial_speed(end-3:end);
                    
                    trial_speed = trial_speed(1:end-4);
                    [~,sidx]=sort(trial_speed,'descend');
                    [~,sidx_gain]=sort(trial_speed_gain,'descend');
                    tmp_rank = 1:18;
                    tmp_rank(sidx)=tmp_rank;
                    tmp_rank_gain = 19:22;
                    tmp_rank_gain(sidx_gain)=tmp_rank_gain;
                    tmp_rank = [tmp_rank, tmp_rank_gain];
                    tmp = data_out{iRep}.allSpikes{iC}(:,2);
                    for iT=1:numel(tmp)
                        tmp(iT)=tmp_rank(tmp(iT));
                    end
                    
                    
                    subplot(2,3,1)
                    scatter(data_out{iRep}.allSpikes{iC}(:,1),tmp,15,get_color(1,contrast),'.')
                    gain_idx = tmp>18;
                    hold on
                    scatter(data_out{iRep}.allSpikes{iC}(gain_idx,1),tmp(gain_idx),15,get_color(gain,contrast),'.')
                    xlim(data_out{iRep}.max_ind(iC)+[-60 40] )
                    ylim([0 23])
                    xline(data_out{iRep}.max_ind(iC))
                    box off
                    title(data_out{iRep}.region{iC})
                    
                    subplot(2,3,2)
                    scatter(trial_speed(sidx),[1:18],15,get_color(1,contrast),'.')
                    
                    hold on
                    scatter(trial_speed_gain(sidx_gain),[19:22],15,get_color(gain,contrast),'.')
                    box off
                    xlabel('cm/s')
                    
                    subplot(2,3,4)
                    scatter(data_out{iRep}.allSpikes{iC}(:,3),tmp,15,get_color(1,contrast),'.')
                    hold on
                    scatter(data_out{iRep}.allSpikes{iC}(gain_idx,3),tmp(gain_idx),15,get_color(gain,contrast),'.')
                    xlabel(sprintf('%.2f',data_out{iRep}.factors(iC)))
                    xlim(data_out{iRep}.max_ind(iC)+[-60 40] )
                    ylim([0 23])
                    xline(data_out{iRep}.max_ind(iC))
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

FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
reg= [];
order = [2,1];
for iF=1:size(output,1)
    for iRep=1:2
    if isempty(output{iF}{iRep})
        disp('help')
        continue
    end
 
    FAST =cat(1,FAST,output{iF}{iRep}.all_fast);
    SLOW = cat(1,SLOW,output{iF}{iRep}.all_slow);
    GAIN = cat(1,GAIN,output{iF}{iRep}.all_gain);
    STAB = cat(1,STAB,output{iF}{iRep}.similarity);
    FACT = cat(1,FACT,output{iF}{order(iRep)}.factors');
    reg = cat(2,reg,output{iF}{iRep}.region);
    end
end


figure
idx = STAB>.4 & ismember(reg,'VISp')';
med_fact = nanmedian(FACT(idx));

idx_1 = idx & FACT<-.1;
idx_2 = idx & FACT>-0.05;
subplot(2,1,1)
plot(x_vec,nanmean(FAST(idx_1,:)))
hold on
plot(x_vec,nanmean(SLOW(idx_1,:)));
plot(x_vec,nanmean(GAIN(idx_1,:)))
legend({'fast','slow','gain'})
grid on
box off
title('delayed cells')
subplot(2,1,2)
plot(x_vec,nanmean(FAST(idx_2,:)))
hold on
plot(x_vec,nanmean(SLOW(idx_2,:)));
plot(x_vec,nanmean(GAIN(idx_2,:)))
legend({'fast','slow','gain'})
grid on
box off
title('non delayed cells')
%% mean across all cells
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
reg= [];
for iF=1:size(output,1)
    
    
    for iRep = 1:2
        if ~isempty(output{iF}{iRep})
            FAST =cat(1,FAST,output{iF}{iRep}.all_fast);
            SLOW = cat(1,SLOW,output{iF}{iRep}.all_slow);
            GAIN = cat(1,GAIN,output{iF}{iRep}.all_gain);
            STAB = cat(1,STAB,output{iF}{iRep}.stability);
            FACT = cat(1,FACT,output{iF}{iRep}.factors');
            reg = cat(2,reg,output{iF}{iRep}.region);
        end
    end
end


figure
region = 'ECT';
idx = STAB>.5 & startsWith(reg,region)';
plot(x_vec,nanmean(FAST(idx,:)))
hold on
plot(x_vec,nanmean(SLOW(idx,:)));
plot(x_vec,nanmean(GAIN(idx,:)))
legend({'fast','slow','gain'})
set(gcf,'Renderer','Painters')
grid on
box off
%saveas(gcf,sprintf('/oak/stanford/groups/giocomo/attialex/FIGURES/peaks_%s.pdf',region))
%% mean across site
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
reg= [];
region = 'VISp';

for iF=1:size(output,1)
    if isempty(output{iF})
        continue
    end
    if numel(output{iF}{1}.stability) ~= numel(output{iF}{1}.region)
        disp(iF)
        continue
    end
    for iRep =1:2
        if ~isempty(output{iF}{iRep})
            idx_region = startsWith(output{iF}{iRep}.region,region)';
            idx = output{iF}{iRep}.similarity>.4 & startsWith(output{iF}{iRep}.region,region)';
            if nnz(idx)>2 && nnz(idx)/nnz(idx_region)>.2
                sel_reg = mean(output{iF}{iRep}.similarity(idx))

                FAST =cat(1,FAST,nanmean(output{iF}{iRep}.all_fast(idx,:),1));
                SLOW = cat(1,SLOW,nanmean(output{iF}{iRep}.all_slow(idx,:),1));
                GAIN = cat(1,GAIN,nanmean(output{iF}{iRep}.all_gain(idx,:),1));
            end
        end
    end
end


figure
boundedline(x_vec,nanmean(FAST),nanstd(FAST)/sqrt(size(FAST,1)),'alpha','cmap',[0 0 0])
hold on
%plot(x_vec,nanmean(FAST))
boundedline(x_vec,nanmean(SLOW),nanstd(FAST)/sqrt(size(FAST,1)),'alpha')
boundedline(x_vec,nanmean(GAIN),nanstd(FAST)/sqrt(size(FAST,1)),'alpha','cmap',get_color(gain,100))


legend({'fast','slow','gain'})
set(gcf,'Renderer','Painters')
grid on
box off

%saveas(gcf,sprintf('/oak/stanford/groups/giocomo/attialex/FIGURES/peaks_%s_new.pdf',region))

%%
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
reg= [];
for iF=1:numel(output)
    if isempty(output{iF})
        continue
    end
    if numel(output{iF}.stability) ~= numel(output{iF}.region)
        disp(iF)
        continue
    end
    FAST =cat(1,FAST,output{iF}.all_fast);
    SLOW = cat(1,SLOW,output{iF}.all_slow);
    GAIN = cat(1,GAIN,output{iF}.all_gain);
    STAB = cat(1,STAB,output{iF}.stability);
    reg = cat(2,reg,output{iF}.subregion);
end


figure
idx = STAB>.2 & ismember(reg,'VISp2/3')';
plot(nanmean(FAST(idx,:)))
hold on
plot(nanmean(SLOW(idx,:)));
plot(nanmean(GAIN(idx,:)))
legend({'fast','slow','gain'})
