%% mean across all cells
FAST =[];
SLOW=[];
GAIN = [];
STAB=[];
FACT = [];
DEPTH = [];
reg= [];
RANGE=[];
SID = [];
cntr = 0;
files = dir('/Volumes/Samsung_T5/attialex/peak_shift22cmFilt_16prec/*.mat');

for iF=1:size(files,1)
    
    
    output = load(fullfile(files(iF).folder,files(iF).name));
    if ~isfield(output,'data_peak')
        continue
    end
    cntr=cntr+1;
            FAST =cat(1,FAST,output.data_peak.all_fast);
            SLOW = cat(1,SLOW,output.data_peak.all_slow);
            GAIN = cat(1,GAIN,output.data_peak.all_gain);
            STAB = cat(1,STAB,output.data_peak.similarity);
            RANGE = cat(1,RANGE,range(output.data_peak.TrialSpeed,2));
            reg = cat(2,reg,output.data_peak.region);
            
            tmp = output.data_peak.depth;
            if isrow(tmp)
                tmp = tmp';
            end
            DEPTH = cat(1,DEPTH,tmp);
            SID = cat(1,SID,ones(size(tmp))*cntr);
    
    %for iRep = 1:2
        
        %if ~isempty(output{iF}{iRep})
%             cntr=cntr+1;
%             FAST =cat(1,FAST,output{iF}{iRep}.all_fast);
%             SLOW = cat(1,SLOW,output{iF}{iRep}.all_slow);
%             GAIN = cat(1,GAIN,output{iF}{iRep}.all_gain);
%             STAB = cat(1,STAB,output{iF}{iRep}.similarity);
%             RANGE = cat(1,RANGE,range(output{iF}{iRep}.TrialSpeed,2));
%             reg = cat(2,reg,output{iF}{iRep}.region);
%             
%             tmp = output{iF}{iRep}.depth;
%             if isrow(tmp)
%                 tmp = tmp';
%             end
%             DEPTH = cat(1,DEPTH,tmp);
%             SID = cat(1,SID,ones(size(tmp))*cntr);
        %end
    %end
end
%%
X=[];
G=[];
x_vec =-20:1:20;
gain = 0.8;
regions = {'MEC','VISp','RS'};
DATA_SLOW = struct();
DATA_FAST = struct();
DATA_GAIN = struct();
N=struct();
figure('Position',[191         677        1049         420])
for iR=1:numel(regions)
    subplot(2,3,iR)
idx = STAB>.5 & startsWith(reg,regions{iR})' & ~isnan(RANGE);
boundedline(x_vec,nanmean(FAST(idx,:)),nanstd(FAST(idx,:))/sqrt(size(FAST(idx,:),1)),'alpha','cmap',[0 0 0])
hold on
%plot(x_vec,nanmean(FAST))
boundedline(x_vec,nanmean(SLOW(idx,:)),nanstd(SLOW(idx,:))/sqrt(size(FAST(idx,:),1)),'alpha')
boundedline(x_vec,nanmean(GAIN(idx,:)),nanstd(GAIN(idx,:))/sqrt(size(FAST(idx,:),1)),'alpha','cmap',get_color(gain,100))
legend({'fast','slow','gain'})
set(gcf,'Renderer','Painters')
grid on
box off
title(regions{iR})
frac = nnz(idx)/nnz(startsWith(reg,regions{iR}));
xlabel(sprintf('%d of %d',nnz(idx),nnz(startsWith(reg,regions{iR}))))
subplot(2,3,4)
hold on
histogram(RANGE(idx),0:5:70,'normalization','probability')
X=cat(1,X,RANGE(idx));
G=cat(1,G,iR*ones(nnz(idx),1));
DATA_SLOW.(regions{iR})=SLOW(idx,:);
DATA_FAST.(regions{iR}) = FAST(idx,:);
DATA_GAIN.(regions{iR}) = GAIN(idx,:);
N.(regions{iR}) = nnz(startsWith(reg,regions{iR}));
end
legend(regions)
%%
V1_speed = sort(X(G==2));
MEC_speed = sort(X(G==1));
RS_speed = sort(X(G==3));
VISp_idx = true(size(V1_speed));
MEC_idx = true(size(MEC_speed));
RS_idx = true(size(RS_speed));

cntr = 0;
while ranksum(V1_speed(VISp_idx),MEC_speed(MEC_idx))<0.05
    cntr = cntr+1;
    VISp_idx(end-cntr+1)=false;
    MEC_idx(cntr)=false;
end
% figure
% histogram(V1_speed(VISp_idx),'normalization','probability')
% hold on
% histogram(MEC_speed(MEC_idx),'normalization','probability')
% figure
% histogram(MEC_speed(MEC_idx))
% hold on
% histogram(MEC_speed)
figure('Position',[ 192         401        1048         269],'Renderer','Painters')

for ii=1:3
    subplot(1,3,ii)
    idx = eval(sprintf('%s_idx',regions{ii}));
    
    boundedline(x_vec,nanmean(DATA_FAST.(regions{ii})(idx,:)),nanstd(DATA_FAST.(regions{ii})(idx,:))/sqrt(nnz(idx)),'alpha','cmap',[0 0 0])
hold on
    boundedline(x_vec,nanmean(DATA_SLOW.(regions{ii})(idx,:)),nanstd(DATA_SLOW.(regions{ii})(idx,:))/sqrt(nnz(idx)),'alpha')
    boundedline(x_vec,nanmean(DATA_GAIN.(regions{ii})(idx,:)),nanstd(DATA_GAIN.(regions{ii})(idx,:))/sqrt(nnz(idx)),'alpha','cmap',get_color(gain,100))
box off
grid on
ylim([0.3 1.15])
title(regions{ii})
xlabel(sprintf('%d of %d',nnz(idx),N.(regions{ii})))
end
    
saveas(gcf,fullfile('/Volumes/Samsung_T5/attialex/images','peakAlignement.pdf'))

