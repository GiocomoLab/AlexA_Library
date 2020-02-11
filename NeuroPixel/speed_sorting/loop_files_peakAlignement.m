
ops.factors = -.25:0.01:.25;
ops.edges = 0:2:400;
ops.search_range=[11:189]; % in bins

ops.midpoints = ops.edges(1:end-1)+(ops.edges(1)+ops.edges(2))/2;
ops.BinWidth = 2;
ops.nBins = numel(ops.edges)-1;

%ops.trials = find(data.trial_gain ==1 & data.trial_contrast==100);
ops.trials = 3:24;
ops.TimeBin = 0.02;
ops.speedWindow = [-10 -1]; % in cm

fi = gausswin(5);
fi=fi'/sum(fi);
ops.filter = fi;
ops.plotfig = false;
ops.idx = [10:2:390]/2;% in bins

OAK='/oak/stanford/groups/giocomo/';

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
[filenames,triggers] = getFilesCriteria('VISp',100,0.5,'F:/NP_DATA');
savepath='F:/temp/images_peak';
if ~isfolder(savepath)
    mkdir(savepath)
end

%%
plotFigures=true;
output = cell(1,numel(filenames));
for iF=1%:numel(filenames)
    try
        %data = load(fullfile(data_dir,filenames{iF}));
        data = load(filenames{iF});
        ops.trials=triggers{iF}(1)+[-18:3];
        data_out = findPeakAlignement(data,ops);
        
        if plotFigures
            x_vec =( -10:10)*ops.BinWidth;
            fig = figure('Position',[680   758   560   220],'visible','off');

            for iC=1:numel(data_out.region)
                if isnan(data_out.max_ind(iC)) || ~(data_out.stability(iC)>.2)
                    continue
                end
                trial_speed = getSpeedAroundPoint(data_out.speed,data.posx,data.trial,ops,data_out.max_ind(iC),ops.speedWindow);
                trial_speed = trial_speed(1:end-4);
                [~,sidx]=sort(trial_speed,'descend');
                tmp_rank = 1:18;
                tmp_rank(sidx)=tmp_rank;
                tmp_rank = [tmp_rank, 19:22];
                tmp = data_out.allSpikes{iC}(:,2);
                for iT=1:numel(tmp)
                    tmp(iT)=tmp_rank(tmp(iT));
                end
                
                
                subplot(2,2,1)
                plot(data_out.allSpikes{iC}(:,1),tmp,'.')
                gain_idx = tmp>18;
                hold on
                plot(data_out.allSpikes{iC}(gain_idx,1),tmp(gain_idx),'r.')
                xlim(data_out.max_ind(iC)+[-60 40] )
                ylim([0 23])
                xline(data_out.max_ind(iC))
                box off
                title(data_out.region{iC})

                subplot(2,2,3)
                plot(data_out.allSpikes{iC}(:,3),tmp,'.')
                hold on
                plot(data_out.allSpikes{iC}(gain_idx,3),tmp(gain_idx),'r.')
                xlabel(sprintf('%.2f',data_out.factors(iC)))
                xlim(data_out.max_ind(iC)+[-60 40] )
                ylim([0 23])
                xline(data_out.max_ind(iC))
                if ~isempty(data_out.subregion)
                    title(data_out.subregion{iC})
                end

                box off
                subplot(2,2,[2 4])
                plot(x_vec,data_out.all_fast(iC,:))
                hold on
                plot(x_vec,data_out.all_slow(iC,:))
                plot(x_vec,data_out.all_gain(iC,:));
                legend({'fast','slow','gain'});
                grid on
                box off
                [~,sn,~]=fileparts(filenames{iF});
                im_save_path = fullfile(savepath,sn);
                if ~isfolder(im_save_path)
                    mkdir(im_save_path)
                end
                impath = fullfile(im_save_path,sprintf('%d.png',data_out.CID(iC)));
                saveas(gcf,impath);
                clf
            end
            close(fig)
        end
        output{iF}=data_out;
    catch ME
        disp(ME.message)

        disp('whoopsie')
    end
end

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
    reg = cat(2,reg,output{iF}.region);
end


figure
idx = STAB>.2 & ismember(reg,'MEC')';
plot(nanmean(FAST(idx,:)))
hold on
plot(nanmean(SLOW(idx,:)));
plot(nanmean(GAIN(idx,:)))
legend({'fast','slow','gain'})
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
idx = STAB>.5 & ismember(reg,'VISp5')';
plot(nanmean(FAST(idx,:)))
hold on
plot(nanmean(SLOW(idx,:)));
plot(nanmean(GAIN(idx,:)))
legend({'fast','slow','gain'})